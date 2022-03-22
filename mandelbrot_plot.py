import matplotlib.pyplot as plt
import pandas as pd

reals = []
imags = []

files = ["build/data_t1.txt", "build/data_t2.txt", "build/data_t3.txt", "build/data_t4.txt"]

for i in files:
   df = pd.read_csv(i, delimiter = "\n", header = None)
   data = list(df[[0]][0])
   data = [x.replace(")", "") for x in data]
   data = [x.replace("(", "") for x in data]

   for j in data:
      x = j.split(",")
      reals.append(float(x[0]))
      imags.append(float(x[1]))


plt.scatter(x=reals, y = imags, marker = ".")
plt.xlabel("Real Axis")
plt.ylabel("Imaginary Axis")
plt.show()