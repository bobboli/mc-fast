import os
import numpy as np

# ./Profile.exe <res> <sphere_radius> <block_X> <block_Y> <block_Z>
main = "Profile.exe"
radius = 0.4
res = [16,32,64]
block_x = 1
block_y = 1

print("runtime in cycles, don't forget to change FREQUENCY in /Profile/Profile.h to your device frequency:")
baseline_result = []
blocking_result = []
vectorization_result = []
for i in res:
    list1 = []
    list2 = []
    list3 = []
    for j in range(30):
        blocking_param = str(block_x)+" "+str(block_y)+" "+str(i)
        f = os.popen(main+" "+str(i)+" "+str(radius)+" "+blocking_param)
        data = f.readlines()
        f.close()
        list1.append(float(data[4].split(' ')[6]))
        list2.append(float(data[6].split(' ')[6]))
        list3.append(float(data[8].split(' ')[6]))
    print("baseline res = "+str(i)+": "+str(np.median(list1)))
    print("blocking res = "+str(i)+": "+str(np.median(list2))+", current blocking by: "+blocking_param)
    print("vectorization res = "+str(i)+": "+str(np.median(list3)))
    baseline_result.append(np.median(list1))
    blocking_result.append(np.median(list2))
    vectorization_result.append(np.median(list3))

print("\nbaseline\n")
for t in baseline_result:
    print(t)

print("\nblocking\n")
for t in blocking_result:
    print(t)

print("\nvectorization\n")
for t in vectorization_result:
    print(t)