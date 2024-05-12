
orthonormal <- function(M) # The set of vectors is orthonormal when the matrix, 
  #which stored these vectors multiplied by the transpose matrix is equal to the identity matrix.
{ 
  if(identical(round(M %*% t(M),8), diag(ncol(M)))) 
   {
    return(TRUE)
   } 
  else 
  {
    return(FALSE)
  }
}

projection <- function(u,v)
{
  proj <- sum(u*v)/(sum(u*u))*u
  return(proj)
}

gramschmidt <- function(M,Mt,i,j,dim) 
{ 
  if(det(M)==0)
  {
    print("The set of vectors is linearly dependent.")
  }else
  {
  Mt[i,]=Mt[i,]-projection(Mt[j,],M[i,])
  if(j<(i-1)) # j is the index of the vectors, which I am projecting on.
  {
    j = j+1
  return (gramschmidt(M,Mt,i,j,dim))
  }
  else if(i<dim) # i is the index of the target vectors.
  { 
    i=i+1 
    return(gramschmidt(M,Mt,i,1,dim))
  }
  else 
    {
      return (round(Mt,8))
    }
  }
}

print("Write dimension: ")
dim <- as.integer(readline()) #I define the dimensionality of space.
M <- matrix(runif(dim * dim,0,10), nrow = dim) # Set of vectors are stored in the matrix. 
#Each row is a different vector.
print("Set of vectors:")
print(M)
print("Are they orthonormal?")
print(orthonormal(M)) 
GS = gramschmidt(M,M,2,1,dim)
print("Matrix after Gram Schmidt orthogonalization")
print(GS)
GS=GS/sqrt(rowSums(GS*GS)) #I normalize the matrix GS by dividing it by the square root of the sum 
#of each squared rows.
print("Normalized matrix:")
print(GS)
print("Is it orthonormal?")
print(orthonormal(GS))