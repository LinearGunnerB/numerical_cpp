#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fname", nargs='+', help="txt file with data in columns like: x y z")
    parser.add_argument("-c","--colors", default='b', nargs='+', help="colors and symbol settings i.e.: -or for dashed red line with o at points")
    parser.add_argument("-g","--grid", default=True, nargs='+', help="turn grid on or off")
    parser.add_argument("-xl","--xlabel", default="", help="xlabel")
    parser.add_argument("-yl","--ylabel", default="", help="ylabel")
    parser.add_argument("-t","--title", default="", help="title")
    parser.add_argument("-l","--legend", nargs='+', default=['y'], help="legend")
    parser.add_argument("-s","--saveas", default="", help="save plot as this file type. Examples\n myplot.svg (for use in html), myplot.eps (for use in pdf)")
    parser.add_argument("-lw","--linewidth", type=int, default=[2], nargs='+', help="")
    parser.add_argument("-ms","--markersize", type=int, default=[2], nargs='+',  help="")
    parser.add_argument("-dpi","--dpi", type=int, default=1600, help="")
    parser.add_argument("-sh","--show", default=True, help="")
    parser.add_argument("-sf","--singlefig", default=True, help="plot each dataset as single figure if True (default), or on the same figure if False")
    args = parser.parse_args()
    
    if len(args.linewidth) > 1 and len(args.markersize) == 1:
        value = args.markersize[0]
        #args.markersize = [value]
        for i in range(len(args.linewidth)-1):
            args.markersize.append(value)

    print(args)
    # Access each variable via args.varname, i.e., args.grid should show True
    return(args)


def plot2D(x,y,args,idx,sf=True):
    h = plt.plot(x,y,args.colors[idx],linewidth=args.linewidth[idx],markersize=args.markersize[idx])
    args.title = ' '.join(args.title)
    plt.title(args.title)
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.grid(args.grid)
    if sf:
        plt.legend(h,args.legend,loc='best')

    return(h)

def main():
    print("\nStarting Python Script")
    args = cmdline()
    print("args.singlefig =",args.singlefig)

    if args.singlefig == True:
        print(">>>> if <<<<<")
        for i in range(len(args.fname)):
            x, y = np.loadtxt(args.fname[i],usecols=(0,1),unpack=True)
            plot2D(x,y,args,i)

            if args.saveas:
                plt.savefig(args.saveas, format=args.saveas.split('.')[1], dpi=args.dpi, bbox_inches='tight')
            if args.show:
                plt.show()
    else:
        print(">>>> else <<<<<<")
        for i in range(len(args.fname)):
            x, y = np.loadtxt(args.fname[i],usecols=(0,1),unpack=True)
            plt.plot(x,y,args.colors[i],linewidth=args.linewidth[i],markersize=args.markersize[i],label=args.legend[i])

        plt.title(args.title)
        plt.xlabel(args.xlabel)
        plt.ylabel(args.ylabel)
        plt.grid(args.grid)
        plt.legend(loc='best')
        
        if args.saveas:
            plt.savefig(args.saveas, format=args.saveas.split('.')[1], dpi=args.dpi, bbox_inches='tight')
        if args.show:
            plt.show()

    print("Ending Python Script")

    
    



if __name__ == "__main__":
    # execute only if run as a script
    main()