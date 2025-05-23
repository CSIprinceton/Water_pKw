# Instructions

1) Compile the code (you need GNU fortran compiler)

```
make
```

2) Make the input file. Create a file (say, input.dat) with the following content

```
natoms nframes  nequil stride nhist
boxa boxb boxc
```

natoms: number of atoms  
nframes: total number of snapshots in your trajectory  
nequil: number of initial equilibration steps to be discarded  
stride: run the analysis at each "stride" steps  
nhist: number of bins used to compute the RDF  
boxa: box length along a  
boxb: box length along b  
boxc: box length along c  

3) Run the code

```
echo a b | ./rdf_wi.x traj.lammpstrj input.dat
```

a and b are 1 or 2, and they correspond to the element index. 1 = O, 2= H. a is the atom belonging to the water ion. b is any other atom.  
traj.lammpstrj is the name of the lammpstrj trajectory file (coordinated have to be in scaled atom coordinates)

Notes: this RDF codes is not general, and is designed to work with the atomic convention adopted in this repo.  
Please modify the source code if you want to make this code applicable to a different systems.
