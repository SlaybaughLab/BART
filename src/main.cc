/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010
 *         Timo Heister, University of Goettingen, 2009, 2010
 */

/* ---------------------------------------------------------------------
 *
 * Author: Weixiong Zheng
 *
 */

#include "../include/model_manager.h"

using namespace dealii;

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;
    
    int dimension;
    if (argc!=3)
    {
      std::cerr << "Call the program as mpirun -np num_proc xtrans input_file_name" << std::endl;
      return 1;
    }
    else
    {
      std::stringstream convert(argv[2]);
      convert >> dimension;
      assert (dimension == 2 || dimension == 3);
    }
    ParameterHandler prm;
    switch (dimension)
    {
      case 2:
      {
        ProblemDefinition<2>::declare_parameters (prm);
        prm.read_input(argv[1]);
        std::string transport_model_name = ProblemDefinition<2>::get_transport_model (prm);
        std::cout << "building done" << std::endl;
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        std_cxx11::shared_ptr<TransportBase<2> > transport_model
        = ModelManager<2>::build_transport_model (transport_model_name, prm);
        //= TransportBase<2>::build_transport_model (transport_model_name,prm);
        transport_model->run ();
        break;
      }
        
      case 3:
      {
        ProblemDefinition<3>::declare_parameters (prm);
        prm.read_input(argv[1]);
        std::string transport_model_name = ProblemDefinition<3>::get_transport_model (prm);
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        std_cxx11::shared_ptr<TransportBase<3> > transport_model
        = ModelManager<3>::build_transport_model (transport_model_name, prm);
        //= TransportBase<3>::build_transport_model (transport_model_name,prm);
        transport_model->run ();
        break;
      }
        
      default:
        AssertThrow (dimension>1,
                     ExcMessage("1D is not implemented for now."));
    }
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
    << "----------------------------------------------------"
    << std::endl;
    std::cerr << "Exception on processing: " << std::endl
    << exc.what() << std::endl
    << "Aborting!" << std::endl
    << "----------------------------------------------------"
    << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
    << "----------------------------------------------------"
    << std::endl;
    std::cerr << "Unknown exception!" << std::endl
    << "Aborting!" << std::endl
    << "----------------------------------------------------"
    << std::endl;
    return 1;
  }
  return 0;
}
