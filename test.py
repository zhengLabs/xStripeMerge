import os

def test():
    for x in range(2,4):
        for m in range(2,4):
            k = 4
            while k <= 64 :
                for n_index in [2,4]:
                    n = n_index * (x * k + m)
                    print(x,"3000",n,k,m)
                    os.system("/home/htf/xRS-merge/build/bin/matching_main %s 3000 %s %s %s" % (x,n,k,m))
                k = k * 2

if __name__ == '__main__':
    test()
