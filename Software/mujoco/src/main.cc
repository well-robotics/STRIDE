//  ros2_mujoco,
// WELL-robotics
// Xiaobin Xiong, modified based on mujoco example and deepbreak example.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "array_safety.h"
#include "glfw_adapter.h"
#include "simulate.h"
#include <mujoco/mujoco.h>

#include "mujoco_sim_helper.hpp"
#include "MuJoCoMessageHandler.h"
// #include "ros_utilities/ros_utilities.hpp"

#define MUJOCO_PLUGIN_DIR "mujoco_plugin"

extern "C"
{
#include <sys/errno.h>
#include <unistd.h>
}

namespace
{
  namespace mj = ::mujoco;
  namespace mju = ::mujoco::sample_util;

  // constants
  const double syncMisalign =
      0.1; // maximum mis-alignment before re-sync (simulation seconds)
  const double simRefreshFraction =
      0.7;                       // fraction of refresh available for simulation
  const int kErrorLength = 1024; // load error string length

  // model and data
  mjModel *m = nullptr;
  mjData *d = nullptr;

  // control noise variables
  mjtNum *ctrlnoise = nullptr;
  std::shared_ptr<mujoco_ros2::MuJoCoMessageHandler::ActuatorCmds> actuator_cmds_ptr;

  using Seconds = std::chrono::duration<double>;

  //---------------------------------------- plugin handling
  //-----------------------------------------

  // return the path to the directory containing the current executable
  // used to determine the location of auto-loaded plugin libraries
  std::string getExecutableDir()
  {

    constexpr char kPathSep = '/';
    const char *path = "/proc/self/exe";
    std::string realpath = [&]() -> std::string
    {
      std::unique_ptr<char[]> realpath(nullptr);
      std::uint32_t buf_size = 128;
      bool success = false;
      while (!success)
      {
        realpath.reset(new (std::nothrow) char[buf_size]);
        if (!realpath)
        {
          std::cerr << "cannot allocate memory to store executable path\n";
          return "";
        }

        std::size_t written = readlink(path, realpath.get(), buf_size);
        if (written < buf_size)
        {
          realpath.get()[written] = '\0';
          success = true;
        }
        else if (written == -1)
        {
          if (errno == EINVAL)
          {
            // path is already not a symlink, just use it
            return path;
          }

          std::cerr << "error while resolving executable path: "
                    << strerror(errno) << '\n';
          return "";
        }
        else
        {
          // realpath is too small, grow and retry
          buf_size *= 2;
        }
      }
      return realpath.get();
    }();

    if (realpath.empty())
    {
      return "";
    }

    for (std::size_t i = realpath.size() - 1; i > 0; --i)
    {
      if (realpath.c_str()[i] == kPathSep)
      {
        return realpath.substr(0, i);
      }
    }

    // don't scan through the entire file system's root
    return "";
  }

  // scan for libraries in the plugin directory to load additional plugins
  void scanPluginLibraries()
  {
    // check and print plugins that are linked directly into the executable
    int nplugin = mjp_pluginCount();
    if (nplugin)
    {
      std::printf("Built-in plugins:\n");
      for (int i = 0; i < nplugin; ++i)
      {
        std::printf("    %s\n", mjp_getPluginAtSlot(i)->name);
      }
    }
    const std::string sep = "/";
    // try to open the ${EXECDIR}/plugin directory
    // ${EXECDIR} is the directory containing the simulate binary itself
    const std::string executable_dir = getExecutableDir();
    
   
    if (executable_dir.empty())
    {
      return;
    }
    const std::string plugin_dir = getExecutableDir() + sep + MUJOCO_PLUGIN_DIR;
    mj_loadAllPluginLibraries(
        plugin_dir.c_str(), +[](const char *filename, int first, int count)
                            {
        std::printf("Plugins registered by library '%s':\n", filename);
        for (int i = first; i < first + count; ++i) {
          std::printf("    %s\n", mjp_getPluginAtSlot(i)->name);
        } });
  }

  //------------------------------------------- simulation
  //-------------------------------------------

  mjModel *LoadModel(const char *file, mj::Simulate &sim)
  {
    // this copy is needed so that the mju::strlen call below compiles
    char filename[mj::Simulate::kMaxFilenameLength];
    mju::strcpy_arr(filename, file);

    RCLCPP_INFO(rclcpp::get_logger("MuJoCo"), "load model from: %s\n", filename);

    // make sure filename is not empty
    if (!filename[0])
    {
      return nullptr;
    }

    // load and compile
    char loadError[kErrorLength] = "";
    mjModel *mnew = 0;
    if (mju::strlen_arr(filename) > 4 &&
        !std::strncmp(filename + mju::strlen_arr(filename) - 4, ".mjb",
                      mju::sizeof_arr(filename) - mju::strlen_arr(filename) +
                          4))
    {
      mnew = mj_loadModel(filename, nullptr);
      if (!mnew)
      {
        mju::strcpy_arr(loadError, "could not load binary model");
      }
    }
    else
    {
      mnew = mj_loadXML(filename, nullptr, loadError,
                        mj::Simulate::kMaxFilenameLength);
      // remove trailing newline character from loadError
      if (loadError[0])
      {
        int error_length = mju::strlen_arr(loadError);
        if (loadError[error_length - 1] == '\n')
        {
          loadError[error_length - 1] = '\0';
        }
      }
    }

    mju::strcpy_arr(sim.loadError, loadError);

    if (!mnew)
    {
      std::printf("%s\n", loadError);
      return nullptr;
    }

    // compiler warning: print and pause
    if (loadError[0])
    {
      // mj_forward() below will print the warning message
      std::printf("Model compiled, but simulation warning (paused):\n  %s\n",
                  loadError);
      sim.run = 0;
    }

    return mnew;
  }

  // PD_CONTROL: to apply control to the motor in mujoco
  void apply_ctrl(mjModel *m, mjData *d)
  {
    // std::cout<< "start applying ctrl"  <<std::endl;
    for (size_t k = 0; k < actuator_cmds_ptr->actuators_name.size(); k++)
    {

      // std::string actuator_joint_name = actuator_cmds_ptr->actuators_name[k];
      // int motor_id = mj_name2id(m, mjOBJ_ACTUATOR,
      //                             actuator_joint_name.c_str());
      // std::string position_name= actuator_joint_name+"_position";
      // std::string velocity_name= actuator_joint_name+"_velocity";
      // int pos_actuator_id = mj_name2id(m, mjOBJ_ACTUATOR,
      //                              position_name.c_str());
      // int vel_actuator_id = mj_name2id(m, mjOBJ_ACTUATOR,
      //                             velocity_name.c_str());
      // if (pos_actuator_id == -1)
      // {
      //   RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),
      //               "not found the name from the received message in mujoco");
      //   continue;
      // }
      // int pos_sensor_id = mj_name2id(m, mjOBJ_SENSOR, (actuator_joint_name + "_pos").c_str());
      // int vel_sensor_id = mj_name2id(m, mjOBJ_SENSOR, (actuator_joint_name + "_vel").c_str());
      // d->ctrl[pos_actuator_id] = actuator_cmds_ptr->pos_des[k];
      //                           // actuator_cmds_ptr->kp[k] *
      //                           //  (actuator_cmds_ptr->pos_des[k] -
      //                           //   d->sensordata[m->sensor_adr[pos_sensor_id]]);
      // d->ctrl[vel_actuator_id] = actuator_cmds_ptr->vel_des[k];
      //                           //   actuator_cmds_ptr->kd[k] *
      //                           //  (actuator_cmds_ptr->vel_des[k] -
      //                           //   d->sensordata[m->sensor_adr[vel_sensor_id]]);
      // d->ctrl[motor_id] =        actuator_cmds_ptr->torque_feedforward[k];


      std::string actuator_joint_name = actuator_cmds_ptr->actuators_name[k];
      int actuator_id = mj_name2id(m, mjOBJ_ACTUATOR,
                                   actuator_joint_name.c_str());

      // std::cout << "name is:" << actuator_cmds_ptr->actuators_name[k] << std::endl;
      if (actuator_id == -1)
      {
        RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),
                    "not found the name from the received message in mujoco");
        continue;
      }
      // assume the actuator cmd->actuator name is joint name
      int pos_sensor_id = mj_name2id(m, mjOBJ_SENSOR, (actuator_joint_name + "_pos").c_str());
      int vel_sensor_id = mj_name2id(m, mjOBJ_SENSOR, (actuator_joint_name + "_vel").c_str());
      // apply feedforward + PD controller; internally implemented
      // similarly to the QDD actuators
      d->ctrl[actuator_id] = actuator_cmds_ptr->kp[k] *
                                 (actuator_cmds_ptr->pos_des[k] -
                                  d->sensordata[m->sensor_adr[pos_sensor_id]]) +
                             actuator_cmds_ptr->kd[k] *
                                 (actuator_cmds_ptr->vel_des[k] -
                                  d->sensordata[m->sensor_adr[vel_sensor_id]]) +
                             actuator_cmds_ptr->torque_feedforward[k];

      // use xml torque bounds
      double torque_min = m->actuator_ctrlrange[actuator_id * 2];
      double torque_max = m->actuator_ctrlrange[actuator_id * 2 + 1];
      //  std::cout << "actuator_name  " << actuator_joint_name
      //   << "  min-max" << torque_min << " " << torque_max <<
      //   " torque " << d->ctrl[actuator_id] <<  std::endl;
      
      d->ctrl[actuator_id] = // apply torque bounds
          std::min(std::max(torque_min, d->ctrl[actuator_id]), torque_max);

      // std::cout << d->ctrl[actuator_id]  << std::endl;
    }
  }

  // simulate in background thread (while rendering in main thread)
  void PhysicsLoop(mj::Simulate &sim)
  {
    // cpu-sim syncronization point
    std::chrono::time_point<mj::Simulate::Clock> syncCPU;
    mjtNum syncSim = 0;
    // run until asked to exit
    while (!sim.exitrequest.load())
    {
      if (!rclcpp::ok())
      {
        sim.exitrequest.store(true);
      }

      if (sim.droploadrequest.load())
      {
        mjModel *mnew = LoadModel(sim.dropfilename, sim);
        sim.droploadrequest.store(false);

        mjData *dnew = nullptr;
        if (mnew)
          dnew = mj_makeData(mnew);
        if (dnew)
        {
          sim.load(sim.dropfilename, mnew, dnew);

          mj_deleteData(d);
          mj_deleteModel(m);

          m = mnew;
          d = dnew;
          mj_forward(m, d);

          // allocate ctrlnoise
          free(ctrlnoise);
          ctrlnoise = (mjtNum *)malloc(sizeof(mjtNum) * m->nu);
          mju_zero(ctrlnoise, m->nu);
        }
      }

      if (sim.uiloadrequest.load())
      {
        sim.uiloadrequest.fetch_sub(1);
        mjModel *mnew = LoadModel(sim.filename, sim);
        mjData *dnew = nullptr;
        if (mnew)
          dnew = mj_makeData(mnew);
        if (dnew)
        {
          sim.load(sim.filename, mnew, dnew);

          mj_deleteData(d);
          mj_deleteModel(m);

          m = mnew;
          d = dnew;
          mj_forward(m, d);

          // allocate ctrlnoise
          free(ctrlnoise);
          ctrlnoise = static_cast<mjtNum *>(malloc(sizeof(mjtNum) * m->nu));
          mju_zero(ctrlnoise, m->nu);
        }
      }

      // sleep for 1 ms or yield, to let main thread run
      //  yield results in busy wait - which has better timing but kills battery
      //  life
      if (sim.run && sim.busywait)
      {
        std::this_thread::yield();
      }
      else
      {
        // std::cout << "sleeping for  millsecond" << std::endl;
        // std::this_thread::sleep_for(std::chrono::microseconds(200));
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

      {
        // lock the sim mutex
        const std::lock_guard<std::mutex> lock(sim.mtx);

        // run only if model is present
        if (m)
        {
          // running
          if (sim.run)
          {
            // record cpu time at start of iteration
            const auto startCPU = mj::Simulate::Clock::now();

            // elapsed CPU and simulation time since last sync
            const auto elapsedCPU = startCPU - syncCPU;
            double elapsedSim = d->time - syncSim;

            // inject noise
            if (sim.ctrlnoisestd)
            {
              // convert rate and scale to discrete time (Ornstein–Uhlenbeck)
              mjtNum rate = mju_exp(-m->opt.timestep /
                                    mju_max(sim.ctrlnoiserate, mjMINVAL));
              mjtNum scale = sim.ctrlnoisestd * mju_sqrt(1 - rate * rate);

              for (int i = 0; i < m->nu; i++)
              {
                // update noise
                ctrlnoise[i] =
                    rate * ctrlnoise[i] + scale * mju_standardNormal(nullptr);

                // apply noise
                d->ctrl[i] += ctrlnoise[i];
              }
            }

            // requested slow-down factor
            double slowdown = 100 / sim.percentRealTime[sim.realTimeIndex];

            // misalignment condition: distance from target sim time is bigger
            // than syncmisalign
            bool misaligned = mju_abs(Seconds(elapsedCPU).count() / slowdown -
                                      elapsedSim) > syncMisalign;

            // out-of-sync (for any reason): reset sync times, step
            if (elapsedSim < 0 || elapsedCPU.count() < 0 ||
                syncCPU.time_since_epoch().count() == 0 || misaligned ||
                sim.speedChanged)
            {
              // re-sync
              syncCPU = startCPU;
              syncSim = d->time;
              sim.speedChanged = false;

              // clear old perturbations, apply new
              mju_zero(d->xfrc_applied, 6 * m->nbody);
              sim.applyposepertubations(0); // move mocap bodies only
              sim.applyforceperturbations();

              apply_ctrl(sim.m, sim.d);

              // custom helping function here:
              mujoco_helper::sim_helper_main(m, d);
              // run single step, let next iteration deal with timing
              mj_step(m, d);
            }

            // in-sync: step until ahead of cpu
            else
            {
              bool measured = false;
              mjtNum prevSim = d->time;

              double refreshTime = simRefreshFraction / sim.refreshRate;

              // step while sim lags behind cpu and within refreshTime
              while (Seconds((d->time - syncSim) * slowdown) <
                         mj::Simulate::Clock::now() - syncCPU &&
                     mj::Simulate::Clock::now() - startCPU <
                         Seconds(refreshTime))
              {
                // measure slowdown before first step
                if (!measured && elapsedSim)
                {
                  sim.measuredSlowdown =
                      std::chrono::duration<double>(elapsedCPU).count() /
                      elapsedSim;
                  measured = true;
                }

                // clear old perturbations, apply new
                mju_zero(d->xfrc_applied, 6 * m->nbody);
                sim.applyposepertubations(0); // move mocap bodies only
                sim.applyforceperturbations();

                apply_ctrl(sim.m, sim.d);
                mujoco_helper::sim_helper_main(m, d);
                // call mj_step
                mj_step(m, d);

                // break if reset
                if (d->time < prevSim)
                {
                  break;
                }
              }
            }

            //
          }

          // paused
          else
          {
            // apply pose perturbation
            sim.applyposepertubations(1); // move mocap and dynamic bodies

            // run mj_forward, to update rendering and joint sliders
            mj_forward(m, d);
          }
        }
      } // release std::lock_guard<std::mutex>
    }
  }
} // namespace

//-------------------------------------- physics_thread
//--------------------------------------------

void PhysicsThread(mj::Simulate *sim, const char *filename)
{
  // request loadmodel if file given (otherwise drag-and-drop)
  if (filename != nullptr)
  { //**** it seems it never went into this block
    m = LoadModel(filename, *sim);
    if (m)
      d = mj_makeData(m);
    if (d)
    {
      sim->load(filename, m, d);
      mj_forward(m, d);

      // allocate ctrlnoise
      free(ctrlnoise);
      ctrlnoise = static_cast<mjtNum *>(malloc(sizeof(mjtNum) * m->nu));
      mju_zero(ctrlnoise, m->nu);
    }
    std::cout << "initializeing 0" << std::endl;
  }
  PhysicsLoop(*sim);

  rclcpp::shutdown();

  // delete everything we allocated
  free(ctrlnoise);
  mj_deleteData(d);
  mj_deleteModel(m);
}

//------------------------------------------ main
// run event loop
int main(int argc, const char **argv)
{

  rclcpp::init(argc, argv);
  // print version, check compatibility
  std::printf("MuJoCo version %s\n", mj_versionString());
  if (mjVERSION_HEADER != mj_version())
  {
    mju_error("Headers and library have different versions");
  }

  // scan for libraries in the plugin directory to load additional plugins
  scanPluginLibraries();

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"), "scanPluginLibraries");

  // simulate object encapsulates the UI
  auto sim =
      std::make_unique<mj::Simulate>(std::make_unique<mj::GlfwAdapter>());

  // ros2 node object
  auto message_handler_node =
      std::make_shared<mujoco_ros2::MuJoCoMessageHandler>(sim.get()); // this is the node

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"), "MuJoCoMessageHandler");

  // set up sim parameters.
  message_handler_node->declare_parameter("T_hold_in_air", mujoco_helper::T_hold_in_air);
  mujoco_helper::T_hold_in_air = message_handler_node->get_parameter("T_hold_in_air").as_double();

  message_handler_node->declare_parameter("T_partial_release", mujoco_helper::T_partial_release);
  mujoco_helper::T_partial_release = message_handler_node->get_parameter("T_partial_release").as_double();

  message_handler_node->declare_parameter("T_release_horizontal", mujoco_helper::T_release_horizontal);
  mujoco_helper::T_release_horizontal = message_handler_node->get_parameter("T_release_horizontal").as_double();
  bool if_record_video = false;
  message_handler_node->declare_parameter("record_video", if_record_video);
  if_record_video = message_handler_node->get_parameter("record_video").as_bool();

  std::string video_name = "";
  message_handler_node->declare_parameter("video_name", video_name);
  video_name = message_handler_node->get_parameter("video_name").as_string();

  int video_height = 1080;
  int video_width = 1920;
  message_handler_node->declare_parameter("video_height", video_height);
  video_height = message_handler_node->get_parameter("video_height").as_int();

  message_handler_node->declare_parameter("video_width", video_width);
  video_width = message_handler_node->get_parameter("video_width").as_int();

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),"Param declared");


  // set up a spin function, which is put on another thread
  auto spin_func = [](std::shared_ptr<mujoco_ros2::MuJoCoMessageHandler> node_ptr)
  {
    rclcpp::spin(node_ptr); // spin constantly, the message frequencey is set based on timmer;
  };

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),"node_ptr created");

  // start physics thread
  std::thread physicsthreadhandle(&PhysicsThread, sim.get(), nullptr);

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),"physicsthreadhandle created");

  actuator_cmds_ptr = message_handler_node->get_actuator_cmds_ptr();

  RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),"message_handler_node created");

  // current logic is: if sim is pause, continue publish
  // sim is running separately.

  auto spin_thread = std::thread{spin_func, message_handler_node};

  // start simulation UI loop (blocking call)
  // multi-thread setup
  sim->renderloop(if_record_video, video_name, video_width, video_height);
  spin_thread.join();
  physicsthreadhandle.join();

  return 0;
}