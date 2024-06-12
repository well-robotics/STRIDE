import matlab.engine
import numpy as np
import time
eng = matlab.engine.start_matlab()
start = time.time() 
a = matlab.double([0,0,0,0,0,0,0])
result = eng.feval("pLeftToe",a)
result2 = eng.feval("pRightToe",a)
end = time.time()
print(end-start)