 #!/usr/bin/env python
import matplotlib.pyplot as plt
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
    plt.title('Core / timing')
    plt.xlabel('cores used')
    plt.ylabel('elapsed time')
    plt.savefig('output/timings.png')
    return

def main():
    df = read_timings_txt()
    df = df.T.reset_index(drop=False).T
    df = df.astype(float)
    draw_graph(df)

if __name__ == "__main__":
    main()
