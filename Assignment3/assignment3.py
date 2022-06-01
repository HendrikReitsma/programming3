 #!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import pandas as pd

def read_timings_txt():
    """
    reads timings
    """
    df = pd.read_csv('output/timings.txt', sep='\t')
    return df


def draw_graph(df):
    """
    draws graph
    """
    plt.plot(df.iloc[:,0], df.iloc[:,1])
    plt.savefig('timings.png')
    return

def main():
    df = read_timings_txt()
    draw_graph(df)

if __name__ == "__main__":
    main()
