import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("build/data_matrix.txt", delimiter = "\n", header = None)
data = list(df[[0]][0])
data = [x.replace(")", "") for x in data]
data = [x.replace("(", "") for x in data]
md = data[1:5]

reals = []
imags = []

for i in data:
   x = i.split(",")
   reals.append(float(x[0]))
   imags.append(float(x[1]))

df = pd.read_csv("build/bool_matrix.txt", delimiter = "\n", header = None)
bools = list(df[[0]][0])

plt.scatter(x=reals, y = imags, alpha=bools, marker = ".")
plt.xlabel("Real Axis")
plt.ylabel("Imaginary Axis")
plt.show()