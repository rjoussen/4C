// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_exporter.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


Core::Communication::Exporter::Exporter(MPI_Comm comm)
    : dummymap_(0, 0, comm),
      frommap_(dummymap_),
      tomap_(dummymap_),
      comm_(comm),
      myrank_(Core::Communication::my_mpi_rank(comm)),
      numproc_(Core::Communication::num_mpi_ranks(comm))
{
}

Core::Communication::Exporter::Exporter(
    const Core::LinAlg::Map& frommap, const Core::LinAlg::Map& tomap, MPI_Comm comm)
    : dummymap_(0, 0, comm),
      frommap_(frommap),
      tomap_(tomap),
      comm_(comm),
      myrank_(Core::Communication::my_mpi_rank(comm)),
      numproc_(Core::Communication::num_mpi_ranks(comm))
{
  construct_exporter();
}

void Core::Communication::Exporter::i_send(const int frompid, const int topid, const char* data,
    const int dsize, const int tag, MPI_Request& request) const
{
  if (my_pid() != frompid) return;
  MPI_Isend((void*)data, dsize, MPI_CHAR, topid, tag, get_comm(), &request);
}

void Core::Communication::Exporter::i_send(const int frompid, const int topid, const int* data,
    const int dsize, const int tag, MPI_Request& request) const
{
  if (my_pid() != frompid) return;
  MPI_Isend((void*)data, dsize, MPI_INT, topid, tag, get_comm(), &request);
}

void Core::Communication::Exporter::i_send(const int frompid, const int topid, const double* data,
    const int dsize, const int tag, MPI_Request& request) const
{
  if (my_pid() != frompid) return;
  MPI_Isend((void*)data, dsize, MPI_DOUBLE, topid, tag, get_comm(), &request);
}

void Core::Communication::Exporter::receive_any(
    int& source, int& tag, std::vector<char>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, get_comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_CHAR, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_CHAR, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::receive(
    const int source, const int tag, std::vector<char>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, get_comm(), &status);
  MPI_Get_count(&status, MPI_CHAR, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_CHAR, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::receive_any(
    int& source, int& tag, std::vector<int>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, get_comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_INT, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_INT, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::receive(
    const int source, const int tag, std::vector<int>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, get_comm(), &status);
  MPI_Get_count(&status, MPI_INT, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_INT, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::receive_any(
    int& source, int& tag, std::vector<double>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, get_comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_DOUBLE, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::receive(
    const int source, const int tag, std::vector<double>& recvbuff, int& length) const
{
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, get_comm(), &status);
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(recvbuff.data(), length, MPI_DOUBLE, source, tag, get_comm(), &status);
}

void Core::Communication::Exporter::allreduce(
    std::vector<int>& sendbuff, std::vector<int>& recvbuff, MPI_Op mpi_op)
{
  int length = (int)sendbuff.size();
  if (length > (int)recvbuff.size()) recvbuff.resize(length);

  MPI_Allreduce(sendbuff.data(), recvbuff.data(), length, MPI_INT, mpi_op, get_comm());
}

void Core::Communication::Exporter::broadcast(
    const int frompid, std::vector<char>& data, const int tag) const
{
  int length = static_cast<int>(data.size());
  MPI_Bcast(&length, 1, MPI_INT, frompid, get_comm());
  if (my_pid() != frompid)
  {
    data.resize(length);
  }
  MPI_Bcast((void*)data.data(), length, MPI_CHAR, frompid, get_comm());
}

void Core::Communication::Exporter::construct_exporter()
{
  if (source_map().same_as(target_map())) return;

  // allocate a sendplan array and init to zero
  // send_plan():
  // send_plan()(lid,proc)    = 1 for data with local id lid needs sending to proc
  // send_plan()(lid,proc)    = 0 otherwise
  // send_plan()(lid,MyPID()) = 0 always! (I never send to myself)
  send_plan().resize(num_proc());

  // To build these plans, everybody has to communicate what he has and wants:
  // bundle this info to save on communication:
  int sizes[2];
  sizes[0] = source_map().num_my_elements();
  sizes[1] = target_map().num_my_elements();
  const int sendsize = sizes[0] + sizes[1];
  std::vector<int> sendbuff;
  sendbuff.reserve(sendsize);
  std::copy(source_map().my_global_elements(),
      source_map().my_global_elements() + source_map().num_my_elements(),
      std::back_inserter(sendbuff));
  std::copy(target_map().my_global_elements(),
      target_map().my_global_elements() + target_map().num_my_elements(),
      std::back_inserter(sendbuff));

  for (int proc = 0; proc < num_proc(); ++proc)
  {
    int recvsizes[2];
    recvsizes[0] = sizes[0];
    recvsizes[1] = sizes[1];
    Core::Communication::broadcast(recvsizes, 2, proc, get_comm());
    const int recvsize = recvsizes[0] + recvsizes[1];
    std::vector<int> recvbuff(recvsize);
    if (proc == my_pid()) std::copy(sendbuff.begin(), sendbuff.end(), recvbuff.data());
    Core::Communication::broadcast(recvbuff.data(), recvsize, proc, get_comm());
    const int* want = recvbuff.data() + recvsizes[0];  // this is what proc needs

    // Loop what proc wants and what I have (send_plan)
    if (proc != my_pid())
    {
      for (int i = 0; i < recvsizes[1]; ++i)
      {
        const int gid = want[i];
        if (source_map().my_gid(gid))
        {
          const int lid = source_map().lid(gid);
          send_plan()[proc].insert(lid);
        }
      }
    }
    Core::Communication::barrier(get_comm());
  }
}

void Core::Communication::Exporter::generic_export(ExporterHelper& helper)
{
  if (send_plan().size() == 0) return;
  // if (SourceMap().SameAs(TargetMap())) return;

  //------------------------------------------------ do the send/recv loop
  for (int i = 0; i < num_proc() - 1; ++i)
  {
    int tproc = my_pid() + 1 + i;
    int sproc = my_pid() - 1 - i;
    if (tproc < 0) tproc += num_proc();
    if (sproc < 0) sproc += num_proc();
    if (tproc > num_proc() - 1) tproc -= num_proc();
    if (sproc > num_proc() - 1) sproc -= num_proc();
    // cout << "Proc " << MyPID() << " tproc " << tproc << " sproc " << sproc << endl;
    // fflush(stdout);

    //------------------------------------------------ do sending to tproc
    // gather all objects to be send
    Core::Communication::PackBuffer sendblock;
    std::vector<int> sendgid;
    sendgid.reserve(send_plan()[tproc].size());

    for (int lid : send_plan()[tproc])
    {
      const int gid = source_map().gid(lid);
      if (helper.pack_object(gid, sendblock)) sendgid.push_back(gid);
    }

    // send tproc no. of chars tproc must receive
    std::vector<int> snmessages(2);
    snmessages[0] = sendblock().size();
    snmessages[1] = sendgid.size();

    MPI_Request sizerequest;
    i_send(my_pid(), tproc, snmessages.data(), 2, 1, sizerequest);

    // do the sending of the objects
    MPI_Request sendrequest;
    i_send(my_pid(), tproc, sendblock().data(), sendblock().size(), 2, sendrequest);

    MPI_Request sendgidrequest;
    i_send(my_pid(), tproc, sendgid.data(), sendgid.size(), 3, sendgidrequest);

    //---------------------------------------- do the receiving from sproc
    // receive how many messages I will receive from sproc
    std::vector<int> rnmessages(2);
    int source = sproc;
    int length = 0;
    int tag = 1;
    // do a blocking specific receive
    receive(source, tag, rnmessages, length);
    if (length != 2 or tag != 1) FOUR_C_THROW("Messages got mixed up");

    // receive the objects
    std::vector<char> recvblock(rnmessages[0]);
    tag = 2;
    receive_any(source, tag, recvblock, length);
    if (tag != 2) FOUR_C_THROW("Messages got mixed up");

    // receive the gids
    std::vector<int> recvgid(rnmessages[1]);
    tag = 3;
    receive_any(source, tag, recvgid, length);
    if (tag != 3) FOUR_C_THROW("Messages got mixed up");

    int j = 0;

    UnpackBuffer buffer(recvblock);
    while (!buffer.at_end())
    {
      int gid = recvgid[j];
      helper.unpack_object(gid, buffer);
      j += 1;
    }

    //----------------------------------- do waiting for messages to tproc to leave
    wait(sizerequest);
    wait(sendrequest);
    wait(sendgidrequest);

    // make sure we do not get mixed up messages as we use wild card receives here
    Core::Communication::barrier(get_comm());
  }

  helper.post_export_cleanup(this);
}

void Core::Communication::Exporter::do_export(
    std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>& data)
{
  AnyObjectExporterHelper<Core::LinAlg::SerialDenseMatrix> helper(data);
  generic_export(helper);
}

FOUR_C_NAMESPACE_CLOSE
