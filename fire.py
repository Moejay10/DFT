#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set()
def compound_Interest(Start, Interest, Contributions, Years):
    """Calculates the compound interest

    Args:
        Start (float): The initial amount
        Contributions (float): Yearly contributions
        Interest (float): Yearly interest 
        Years (float): Number of years 
    """
    years = np.linspace(0, Years, Years)
    R = np.linspace(0, Interest, int(Interest*100)+1)
    x = []
    F = {}
    for r in R:
        F[str(int(r*100)) + "%"] = [Start]
    
    for r in R:
        i = 0
        for y in range(len(years)):
            temp = F[str(int(r*100)) + "%"][i]*(1+r) + Contributions
            F[str(int(r*100)) + "%"].append(temp)
            i += 1


        x.append(str(int(r*100)) + "%")
    
    df = pd.DataFrame(F)
    target = np.ones(len(years))*1E7
    df.plot(style='.-')
    plt.plot(years, target, '--k')
    plt.xlabel("Years")
    plt.ylabel("Money (kr)")
    plt.title("Financial Independence Retire Early")
    plt.show() 


if __name__ == '__main__':
    compound_Interest(2E5, 0.1, 2.5E5, 20)

