Saves array A as binary file with name Aname.

````julia
function save_array(Aname,A)
    fname = string(Aname,".bin")
    out = open(fname,"w"); write(out,A); close(out)
end
````

````
save_array (generic function with 1 method)
````

    load_array(Aname,A)

Loads bin files with name Aname

````julia
function load_array(Aname,A)
    fname = string(Aname,".bin")
    fid=open(fname,"r"); read!(fid,A); close(fid)
end
````

````
load_array (generic function with 1 method)
````

    main function

````julia
function main()
    # parameter
    n = 3
    A = rand(Float64,n,n)
    B = zeros(Float64,n,n)

    # save array A
    save_array("../test/LitTest",A)

    # load array A into B
    load_array("../test/LitTest",B)

    return B
end
````

````
main (generic function with 1 method)
````

B = main()
heatmap(B)

