import numpy as np
import matplotlib.pyplot as plt

def read_data(file):
    time = list()
    x = list()
    y = list()
    z = list()

    f = open (file, "r")
    for row in f:
        try:
            data = row.split(',')
            time.append(float(data[0]))
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
        except:
            pass
    return time, x, y, z

time, x, y, z = read_data("yxyx/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 18_41_00.csv")
#time_cor, x_c, y_c, z_c = read_data("SunSensorBody-Combined-Short-2024-07-08 14_15_16.csv_corrected.csv")



plt.plot(time, x, label="original x")
#plt.plot(time_cor, x_c, label="corrected x")
plt.plot(time, y, label="original y")
#plt.plot(time_cor, y_c, label="corrected y")
plt.plot(time, z, label="original z")
#plt.plot(time_cor, z_c, label="corrected z")
plt.xlabel("time")
plt.ylabel("x")
plt.legend()
plt.show()
