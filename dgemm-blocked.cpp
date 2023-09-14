const char* dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   // insert your code here
   int N = n/block_size;
   int b = block_size * block_size;
   for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         double tempc[b];
         for(int u=0;u<block_size;u++){
            for(int h=0;h<block_size;h++){
               tempc[u+h*block_size] = C[i*block_size+u+(j*block_size+h)*n];
            }
         }
         for(int k=0;k<N;k++){
            double tempa[b], tempb[b];
            for(int x=0;x<block_size;x++){
               for(int y=0;y<block_size;y++){
                  tempa[x+y*block_size] = A[x+i*block_size+(k*block_size+y)*n];
                  tempb[x+y*block_size] = B[x+k*block_size+(j*block_size+y)*n];
               }
            }
            for(int l=0;l<block_size;l++){
               for(int m=0;m<block_size;m++){
                  for(int p=0; p<block_size;p++)
                  tempc[l+m*block_size] += tempa[l+p*block_size] * tempb[p+m*block_size];
               }
            } 
            
         }
         for(int u=0;u<block_size;u++){
            for(int h=0;h<block_size;h++){
               C[i*block_size+u+(j*block_size+h)*n] = tempc[u+h*block_size];
            }
         }
      }
   }
}