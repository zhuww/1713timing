import numpy as np
import matplotlib.pyplot as plt
import triangle

data = np.loadtxt('chain_1.txt')

print data.shape

data = data[:,:-4]
figure = triangle.corner(data,  
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True, title_args={'fontsize': 12})
#figure.gca().annotate("A title", xy=(0.5, 1.0), xycoords='figure fraction', 
        #xytext=(0, -5), textcoords="offset points",
        #ha="center", va="top")
figure.savefig('triangle.png')
