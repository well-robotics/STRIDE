template <class T>
class MovingAverageFilter{
  // this filter outputs the averaged reading on window n for data type T
  // requisite: the + and / operators between themselves and with scalar values should be defined
public: 
  int window_size; 
  T previous_data[50];
  T zero;  
  
  MovingAverageFilter(int size, T init){
    this->window_size = size;
    this->zero = init;
    // initialize the array
    int i; 
    for (i = 0; i<this->window_size; i++){ 
      this->previous_data[i] = init; 
    }
  }

  T calc(){
    T output = this->zero; 
    for(int i=0;i<window_size;i++){
      output = output + previous_data[i]; 
    }
    return (output/window_size); 
  }
  void renew(T input){
    for(int i=0;i<window_size;i++){
      previous_data[i]= previous_data[i+1]; 
    }
    previous_data[window_size - 1] = input; 
  }

};
















// class MovingAverageFilter{
// // this filter outputs the averaged reading on window n    
// public:  
//   int window_size; 
//   double previous_data[50] = {0}; // the max size of the window is set to 50
//   MovingAverageFilter(int size){
//     this->window_size = size; 
//   } 
//   double calc(){
//     double output =0; 
//     for(int i=0;i<window_size;i++){
//       output = output + previous_data[i]; 
//     }
//     return (output/window_size); 
//   }
//   void renew(double input){
//     for(int i=0;i<window_size;i++){
//       previous_data[i]= previous_data[i+1]; 
//     }
//     previous_data[window_size - 1] = input; 
//   }
// };