//  #include <iostream>
// #include <stdio.h> // for gcc >= 4.4 compatibility
// #include "CommonMPI.hpp"

// namespace MPI 
// {
// // give single processor initial values to these
// // in practice, the user should call initMpiStuff
// // near the start of there code, after MPI_Init(*,*)...
// int mpi_rank;
// int mpi_size = 1;
// int row_rank;
// int row_size;
// int ndims = dim; // Number of dimensions in domain decomposition
// int reorder = 1;
// int n_rhs = 0;
// int n_abv = 0;
// int colour;
// int *mpi_coords;
// int *n_coords;
// int *n_rank;
// int *periodic;
// int *mpi_dims;

// MPI_Comm mpi_comm = MPI_COMM_WORLD;
// MPI_Comm row_comm;
// MPI_Comm cartcomm;

// MPI_Request request;
// MPI_Request request2;

// MPI_Status status;

// 	void initMpiStuff() 
// 	{	MPI_Comm_rank(mpi_comm, &mpi_rank);
// 	  	MPI_Comm_size(mpi_comm, &mpi_size);
// 	  	MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
// 	}

// 	void initMpiDomain()
// 	{	mpi_dims = new int[2];
// 		mpi_coords = new int[2];
// 		periodic = new int[2];

// 		// now we need to define periodic depending on input file
// 		if (rightBnd.type == 1 && leftBnd.type == 1)
// 		{ 	periodic[0] = 1;
// 		}
// 		if (rightBnd.type != 1 || leftBnd.type != 1)
// 		{ 	periodic[0] = 0;
// 		}
// 		if (topBnd.type == 1 && bottomBnd.type == 1)
// 		{ 	periodic[1] = 1;
// 		}	
// 		if (topBnd.type != 1 || bottomBnd.type != 1)
// 		{ 	periodic[1] = 0;
// 		}
// 		MPI_Dims_create(mpi_size, ndims, mpi_dims);
// 		MPI_Cart_create(mpi_comm, 2, mpi_dims, periodic, reorder, &cartcomm);
// 		// MPI_Cart_coords(cartcomm, mpi_rank, ndims, mpi_coords);
//   	}

//   	// Split Communicators 0 is x, 1 is y
//   	void splitComm(int i)
//   	{	colour = mpi_rank / mpi_dims[i];

// 		// Split the communicator based on the color and use the
// 		// original rank for ordering
		
// 		MPI_Comm_split(MPI_COMM_WORLD, colour, mpi_rank, &row_comm);

		
// 		MPI_Comm_rank(row_comm, &row_rank);
// 		MPI_Comm_size(row_comm, &row_size);

//   	}

// 	void exchangeData(double **dat, int d1, int d2)
// 	{	getNeighbours();
// 		int ng_rhs[mpi_size]; int ngi_rhs[mpi_size];
// 		// int *ng_rhs = new int[mpi_size];
// 	    MPI_Request request;
// 	    MPI_Status status;
// 		int buffer[d2]; 
// 		int buffer2[d2]; 
// 		int left, right;
// 		// cout << n_rhs << "  " << mpi_rank << endl;
// 		for (int i = 0; i < mpi_size; ++i)
// 		{	ng_rhs[i] = 0;
// 			ngi_rhs[i] = i;
// 		}
// 		ng_rhs[mpi_rank] = n_rhs;
// 		splitComm(0);
// 		// cout << mpi_rank << "   " << n_rhs << "   " << n_abv << endl;
// 		MPI_Allgather(&n_rhs, 1, MPI_INT, &ng_rhs, 1, MPI_INT, cartcomm);

// 		// Find the neighbours...
// 	    right = (row_rank + 1) % row_size;
// 	    left = row_rank - 1;
// 	    if (left < 0)
// 	    {	left =  row_size - 1;
// 	    }
// 		// Now fill ghost cells between set up the receives...
// 		// First send (Nx_-1) ---> 0
// 		for (int i = 0; i < d2; ++i)
// 		{	buffer[i] = row_rank;
// 			buffer2[i] = mpi_rank;
// 		}

// 		cout << row_rank << "  " << left << "  " << right << "  " << colour << endl;
// 		cout << row_rank << "  " << buffer[0] << "  " << buffer2[0] << "  " << colour << endl;

// 		// if (mpi_rank==0)
// 		// {	cout << mpi_rank << "  " << buffer[0] << "   2" << endl;
// 		// }
// 		// cout << mpi_rank << "  " << buffer[0] << endl;
// 		MPI_Sendrecv(buffer, d2, MPI_INT, right, 123, buffer2, d2, MPI_INT, left, 123, row_comm, &status);
// 		cout << row_rank << "  " << buffer[0] << "  " << buffer2[0] << "  " << colour << endl;
// 		//First send 1 ---> Nx_
// 		// if (mpi_rank==0)
// 		// {	cout << buffer[0] << "  " << buffer[0] << endl;
// 		// }
// 		// MPI_Irecv(buffer, d2, MPI_Datatype MPI_INT, ng_rhs[mpi_rank], 123, cartcomm, &request);
// 		// MPI_Isend(buffer, d2, MPI_Datatype MPI_INT, ng_rhs[mpi_rank], 123, cartcomm, &request);
// 		// MPI_Waitall(1, &request, &status);

// 		// Now start send
// 		// cout << mpi_rank << "  " << buffer2[0] << endl;






// 		// if (mpi_dims[1]==1)
// 		// Shift array elements mpi_size/2 across so that proc is never sending or receiving itself...
// 		// {	int tmp;	int tmp_i;
// 		// 	for (int it = 0; it < 4; ++it)
// 		// 	{	tmp = ng_rhs[mpi_size-1];
// 		// 		tmp_i = ngi_rhs[mpi_size-1];
// 		// 		for (int i = mpi_size; i != 0; --i)
// 		// 		{	ng_rhs[i] = ng_rhs[i-1];
// 		// 			ngi_rhs[i] = ngi_rhs[i-1];
// 		// 		}
// 		// 		ng_rhs[0] = tmp;
// 		// 		ngi_rhs[0] = tmp_i;
// 		// 	}
// 		// }


// 		// if (mpi_dims[1]==1)
// 		// {	MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[mpi_rank], 123, buffer, d2, MPI_INT, ngi_rhs[mpi_rank], 123, MPI_COMM_WORLD, &status);
// 		// 	// if (mpi_rank < (mpi_size-3))
// 		// 	// {	//MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[mpi_rank+2], 123, buffer, d2, MPI_INT, mpi_rank+2, 123, MPI_COMM_WORLD, &status);
// 		// 	// }
// 		// 	// if (mpi_rank == (mpi_size-3))
// 		// 	// {	//MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[mpi_size-1], 123, buffer, d2, MPI_INT, mpi_size-1, 123, MPI_COMM_WORLD, &status);
// 		// 	// 	cout << mpi_rank << "     " << ng_rhs[mpi_size-1] << "  " << mpi_size-1 << endl;
// 		// 	// }
// 		// 	// if (mpi_rank > (mpi_size-3))
// 		// 	// {	//MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[(mpi_size-1)-mpi_rank], 123, buffer, d2, MPI_INT, (mpi_size-1)-mpi_rank, 123, MPI_COMM_WORLD, &status);
// 		// 	// 	cout << mpi_rank << "     " << ng_rhs[(mpi_size-1)-mpi_rank] << "  " << (mpi_size-1)-mpi_rank << endl;
// 		// 	// }
// 		// }
// 		// if (mpi_dims[1] > 1)
// 		// {	//if (mpi_rank > 0)
// 		// 	// {	MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[mpi_rank-1], 123, buffer, d2, MPI_INT, mpi_rank-1, 123, MPI_COMM_WORLD, &status);
// 		// 	// 	cout << mpi_rank << "     " << ng_rhs[mpi_rank-1] << "  " << mpi_rank-1 << endl;
// 		// 	// }
// 		// 	// else
// 		// 	// {	MPI_Sendrecv(buffer, d2, MPI_INT, ng_rhs[mpi_size-1], 123, buffer, d2, MPI_INT, mpi_size-1, 123, MPI_COMM_WORLD, &status);
// 		// 	// 	cout << mpi_rank << "     " << ng_rhs[mpi_size-1] << "  " << mpi_size-1 << endl;
// 		// 	// }
// 		// }
// 		// cout << buffer[0] << "  " << mpi_rank << endl;
				

// 		// MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_INT, &x, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, mpi_rank, int recvtag,      cartcomm, MPI_Status *status)
// 		// // Find neighbours to the RHS diagonal...
// 		// if (mpi_coords[0] < mpi_dims[0])
// 		// {	rhs[0] = mpi_coords[0]+1;
// 		// 	rhs[1] = mpi_coords[1]+1;
// 		// 	MPI_Cart_rank(cartcomm, rhs, &x);
// 		// 	cout << mpi_rank << "   " << x << endl;

// 		// }
// 		// // Find neighbours to the LHS diagonal...
// 		// if (mpi_coords[0] > 0)
// 		// {	rhs[0] = mpi_coords[0]-1;
// 		// 	rhs[1] = mpi_coords[1]+1;
// 		// 	MPI_Cart_rank(cartcomm, rhs, &x);
// 		// 	cout << mpi_rank << "   " << x << endl;

// 		// }

// 		// MPI_Cart_shift(mpi_comm, 1, 1, 0, &rank);
// 	}

// 	void getNeighbours()
// 	{	// rhs is right hand neighbour, abv is neighbour above current cell...
// 		int *rhs = new int[2];
// 		int *abv = new int[2];
// 		rhs[0] = mpi_coords[0];
// 		rhs[1] = mpi_coords[1];
// 		abv[0] = mpi_coords[0];
// 		abv[1] = mpi_coords[1];

// 		// Find neighbours to the right...
// 		if (mpi_coords[0] < mpi_dims[0])
// 		{	if (mpi_coords[0] == (mpi_dims[0]-1))
// 			{	rhs[0] = 0;
// 				rhs[1] = mpi_coords[1];
// 				MPI_Cart_rank(cartcomm, rhs, &n_rhs);
// 			}
// 			else
// 			{	rhs[0] = mpi_coords[0]+1;
// 				rhs[1] = mpi_coords[1];
// 				MPI_Cart_rank(cartcomm, rhs, &n_rhs);
// 			}
// 		}
// 		// Find neighbours to the above...
// 		if (mpi_coords[1] < mpi_dims[1])
// 		{	if (mpi_coords[1] == (mpi_dims[1]-1))
// 			{	abv[0] = mpi_coords[0];
// 				abv[1] = 0;
// 				MPI_Cart_rank(cartcomm, abv, &n_abv);
// 			}
// 			else
// 			{	abv[0] = mpi_coords[0];
// 				abv[1] = mpi_coords[1]+1;
// 				MPI_Cart_rank(cartcomm, abv, &n_abv);
// 			}
// 		}
// 	}

//   	// Give in the mpi_rank return the domain decomposition co-ordinates
//   	void getRankCoords() 
// 	{	MPI_Cart_coords(cartcomm, mpi_rank, ndims, mpi_coords);
// 	}


// 	// This gets the bounds of the local domains for global referencing
// 	// Bounds (l_bnd) are referenced as follows:
// 	//	 _____[1][1]_____
// 	//  |                |
// 	//  |                |
// 	// [0][0]         [0][1]
// 	//  |                |
// 	//  |                |
// 	//  |_____[1][0]_____|
// 	//
//   	void getGlobalBounds(int **x_bnd, int **y_bnd, int *Nx, int *Ny) 
// 	{ 	int x_rnk[mpi_dims[0]];
// 		int y_rnk[mpi_dims[1]];
// 		// Get x axis order
// 		int putin[2], X;
// 		for (int j = 0; j < mpi_dims[0]; ++j)
// 		{   putin[0] = j;
// 		    putin[1] = 0;
// 		    MPI_Cart_rank(cartcomm, putin, &X);
// 		    // if (mpi_rank == 0)
// 		    {  	x_rnk[j] = X;
// 		    	// cout << j << " " << x_rnk[j] << endl;
// 		    }
// 		}
// 		// sleep(5);
// 		// Get y axis order
// 		for (int j = 0; j < mpi_dims[1]; ++j)
// 		{   putin[0] = 0;
// 		    putin[1] = j;
// 		    MPI_Cart_rank(cartcomm, putin, &X);
// 		    // if (mpi_rank == 0)
// 		    {  y_rnk[j] = X;
// 		    }
// 		}
// 		// Create local x bounds
// 		for (int i = 0; i < mpi_dims[0]; ++i)
// 		{	//if (mpi_rank == 0)
// 			{	if (i == 0)
// 				{	x_bnd[i][0] = 0;
// 					x_bnd[i][1] = Nx[x_rnk[i]];
// 				}
// 				if (i != 0)
// 				{	x_bnd[i][0] = x_bnd[i-1][1] - 2;
// 					x_bnd[i][1] = x_bnd[i][0] + Nx[x_rnk[i]];
// 				}
// 				// x_bnd[i][1] = x_bnd[i-1][1] -2;
// 				// x_bnd[i][1] = x_bnd[i][0];
// 				// x_bnd[i][1] = x_bnd[i][0] + Nx[x_rnk[i-1]];
// 				if (i == (mpi_dims[0]-1))
// 				{	x_bnd[i][1] = Nx_g;
// 				}
// 				// cout << mpi_rank << "  " << x_bnd[i][0] << " " << x_bnd[i][1] << endl;
// 				// cout << i << " " << x_rnk[i] << endl;
// 			}
// 		}
// 		// // Send x bounds to all processors...
// 		// for (int i = 0; i < mpi_dims[0]; ++i)
// 		// {	MPI_Bcast(&x_bnd[i][0], mpi_dims[0], MPI_INT, 0, cartcomm);
// 		// 	MPI_Bcast(&x_bnd[i][1], mpi_dims[0], MPI_INT, 0, cartcomm);
// 		// }
// 		// Create local y bounds
// 		for (int i = 0; i < mpi_dims[1]; ++i)
// 		{	//f (mpi_rank == 0)
// 			{	if (i == 0)
// 				{	y_bnd[i][0] = 0;
// 					y_bnd[i][1] = Ny[y_rnk[i]];
// 				}
// 				if (i != 0)
// 				{	y_bnd[i][0] = y_bnd[i-1][1] - 2;
// 					y_bnd[i][1] = y_bnd[i][0] + Ny[y_rnk[i]];
// 				}
// 				// x_bnd[i][1] = x_bnd[i-1][1] -2;
// 				// x_bnd[i][1] = x_bnd[i][0];
// 				// x_bnd[i][1] = x_bnd[i][0] + Nx[x_rnk[i-1]];
// 				if (i == (mpi_dims[1]-1))
// 				{	y_bnd[i][1] = Ny_g;
// 				}
// 				// cout << mpi_rank << "  " << y_bnd[i][0] << " " << y_bnd[i][1] << endl;
// 				// cout << i << " " << x_rnk[i] << endl;
// 			}
// 		}

// 		// // Send y bounds to all processors...
// 		// for (int i = 0; i < mpi_dims[1]; ++i)
// 		// {	MPI_Bcast(&y_bnd[i][0], mpi_dims[1], MPI_INT, 0, cartcomm);
// 		// 	MPI_Bcast(&y_bnd[i][1], mpi_dims[1], MPI_INT, 0, cartcomm);
// 		// }

// 	}
// 	void getLocalBounds(int **l_bnd, int **x_bnd, int **y_bnd)
// 	{	// Now give each processor his local bounds so we can use it in equations....
// 		for (int i = 0; i < mpi_dims[0]; ++i)
// 		{	for (int j = 0; j < mpi_dims[1]; ++j)
// 			{	if (mpi_coords[0] == i && mpi_coords[1] == j)
// 				{	l_bnd[0][0]	= x_bnd[i][0];
// 					l_bnd[0][1] = x_bnd[i][1];
// 					l_bnd[1][0] = y_bnd[j][0];
// 					l_bnd[1][1] = y_bnd[j][1];
// 				}
// 			}
// 		}
		
// 	}

// }