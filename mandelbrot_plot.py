import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("build/data_matrix.txt", delimiter = "\n", header = None)
data = list(df[[0]][0])
data = [x.replace(")", "") for x in data]
data = [x.replace("(", "") for x in data]


df = pd.read_csv("build/bool_matrix.txt", delimiter = "\n", header = None)
bools = list(df[[0]][0])

reals = []
imags = []

for i in data:
   x = i.split(",")
   reals.append(float(x[0]))
   imags.append(float(x[1]))

reals = [reals[i]*bools[i] for i in range(len(bools))]
imags = [imags[i]*bools[i] for i in range(len(bools))]

plt.scatter(x=reals, y = imags, marker = ".")
plt.xlabel("Real Axis")
plt.ylabel("Imaginary Axis")
plt.show()