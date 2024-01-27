function myfunct(n::Cint, x_::Ptr{Cdouble})
	##########################################
	# This is the Shekel function with m=5
	##########################################
	x = Array(Cdouble,n)
	x = pointer_to_array(x_,n)

	A = Array(Cdouble,5,4)
	for i=1:4
		A[1,i]=4.0
		A[2,i]=1.0
		A[3,i]=8.0
		A[4,i]=6.0
	end
	for i=1:2
		A[5,2*(i-1)+1]=3.0
		A[5,2*i]=7.0
	end

	C = Array(Cdouble,5)
	C[1]=0.10
	C[2]=0.20
	C[3]=0.20
	C[4]=0.40
	C[5]=0.40

	F  = 0.0
	FA = 0.0

	for I=1:5
		for J=1:4
			FA = FA +(x[J]-A[I,J])^2
		end
		if ((FA+C[I])==0.0) then
			F=1.e25
			return convert(Cdouble, F)::Cdouble
		end
		F = F - 1.0/(FA+C[I]) 
		FA = 0.0
	end

	return convert(Cdouble, F)::Cdouble
end
const myfunct_c = cfunction(myfunct, Cdouble, (Cint,Ptr{Cdouble}))

n       = 4
maxiter = 1000
tol     = 1.e-6
lb      = Array(Cdouble,n)
ub      = Array(Cdouble,n)
xin     = Array(Cdouble,n)
xbest   = Array(Cdouble,n)
fbest   = Array(Cdouble,1)

fill!(lb,-10.0)
fill!(ub, 10.0)
fill!(xin, 0.0)
fill!(fbest, 1.0e25)

println("\n\nNow calling SDBOX using funct in julia...")
ccall((:sdbox,"libsdbox.a"),
      Void,
      (Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Void},Ptr{Float64},Ptr{Float64},Int32,Float64),
      n,pointer(lb),pointer(ub),pointer(xin),myfunct_c,pointer(xbest),pointer(fbest),maxiter,tol)

println("\n\nOptimal function value is: ",fbest[1])
println("         Optimal point is: ")
for i = 1:n
	@printf("        xbest[%d] = %13.6e\n",i,xbest[i])
end
