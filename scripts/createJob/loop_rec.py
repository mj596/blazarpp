# Make loop over data in variable length and variable nested level lists 
# How to use it:
# a = your data tuple
# ret = results vector
# Example:
# x=((1,2,3),(4,5,6),(7,8,9))
# iterations=[]
# loop_rec(x,iterations)
# print iterations
# [[7, 4, 1], [7, 4, 2], [7, 4, 3], [7, 5, 1], [7, 5, 2], [7, 5, 3], [7, 6, 1], [7, 6, 2], [7, 6, 3], [8, 4, 1], [8, 4, 2], [8, 4, 3], [8, 5, 1], [8, 5, 2], [8, 5, 3], [8, 6, 1], [8, 6, 2], [8, 6, 3], [9, 4, 1], [9, 4, 2], [9, 4, 3], [9, 5, 1], [9, 5, 2], [9, 5, 3], [9, 6, 1], [9, 6, 2], [9, 6, 3]]
def loop_rec(a,ret,b=None):
    # first make an array consisting of all data in a except for the last one;
    # this will be passed futher
    vpass=[]
    for j in range(len(a)-1):
        vpass.append(a[j])

    # if a has more than one element make a loop over elements in last element
    # which is not present in vpass anymore
    if len(a)>=1:
        for i in a[len(a)-1]:
            # now take all the previous passed data and and current
            argpass=[]
            if b != None:
                for k in b:
                    argpass.append(k)
            argpass.append(i)
            
            # now we can recursively call loop_rec again with vpass and all the history
            # of previous steps
            loop_rec(vpass,ret,argpass)
    else:
        # we are at the end; we can print the data and append it to ret vector
        #        print b
        ret.append(b)
        
