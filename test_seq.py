#! /usr/bin/env python



def run_seq(N=100):
    a  = 3
    b  = 7
    P  = 11
    x0 = 5

    x_cur = x_prev = x0
    for i in range(N):
        print("{:4}: {:3}".format(i, x_cur))
        x_cur = (a*x_prev + b) % P
        x_prev = x_cur
if __name__ == "__main__":
    run_seq()