#declare M=20;
#declare N=20;
#declare Np=M*N;
#declare spins=array[Np];
#declare cam=<M/2,-3*M/2,2*M/3>;
//#declare cam=<M/2,N/2,1.1*M>;
#declare midspace=<M/2,N/2,0>;
#declare i=0;
#while(i<Np)
    #declare spins[i]=0;
    #declare i=i+1;
#end
//-------------
background
{
    color rgb<0,0,0>
}
camera
{
    location cam
    sky<0,0,1>
    right<-1,0,0>
    look_at midspace
}
light_source
{
    cam
    color rgb<1,1,1>
    spotlight radius 1
    adaptive 1
    jitter
    point_at midspace
    shadowless
}
//--------------
#declare up_top=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,-0.5>,<0.5,0.5,0>}
    texture { pigment { color rgb <1,0,0> } }
    finish  { phong 1.0 phong_size 60}
}
#declare up_bottom=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,0>,<0.5,0.5,0.5>}
    texture { pigment { color rgb <1,1,1> } }
    finish  { phong 1.0 phong_size 60}
}
#declare up_spin=
union
{
    object{up_top}
    object{up_bottom}
}
//----------------
#declare down_top=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,-0.5>,<0.5,0.5,0>}
    texture { pigment { color rgb <1,1,1> } }
    finish  { phong 1.0 phong_size 60 }
}
#declare down_bottom=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,0>,<0.5,0.5,0.5>}
    texture { pigment { color rgb <1,0,0> } }
    finish  { phong 1.0 phong_size 60}
}
#declare down_spin=
union
{
    object{down_top}
    object{down_bottom}
}
//----------------

#fopen MyFile "./frame_data.dat" read
#declare i=0;
#while (defined(MyFile))
    #read (MyFile,var1)
    #declare spins[i]=var1;
    #declare i=i+1;
#end
#fclose MyFile

#declare i=0;
#while(i<N)
    #declare j=0;
    #declare jj=i*M;
    #while(j<M)
        #declare ii=jj+j; 
        #if(spins[ii]=1)
            object{up_spin translate<j,i,0>}
        #else
            #if(spins[ii]=-1)
                object{down_spin translate<j,i,0>}
            #end
        #end
    #declare j=j+1;
    #end
#declare i=i+1;
#end


