/*
 * Graph.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <Graph.h>

using namespace std;

namespace sparse {

Graph::Graph() {
  // TODO Auto-generated constructor stub

}

Graph::~Graph() {
  // TODO Auto-generated destructor stub
}

void Graph::addVertex(const Vertex& v_){
  _vertices.push_back(v_);
}

void Graph::addEdge(const Edge& e_){
  _edges.push_back(e_);
}


void Graph::updateVertices(const std::vector<Vertex>& new_vertices_){
  for (int i = 0; i < _vertices.size(); ++i) {
    Vertex new_vertex = new_vertices_[i];
    _vertices[i].setData(new_vertex.data());
  }
}



void Graph::exportToG2OFile(const std::string& filename_) {
  cerr << BOLDYELLOW << "\t" << "Exporting updated graph to: " << filename_ << RESET << endl;
  ofstream file(filename_);

  //! Writing SE3 Vertices
  for (int i = 0; i < _vertices.size(); ++i) {
    Vertex& v = _vertices[i];

    file << VERTEX_SE3 << " " << v.id() << " ";

    Vector3 t = v.data().translation();
    file << t.x() << " " << t.y() << " " << t.z() << " ";

    QuaternionReal q(v.data().linear());
    file << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " " << endl;
  }

  file << endl << endl << endl;

  //! Writing SE3 Edge
  for (int i = 0; i < _edges.size(); ++i) {
    Edge& e = _edges[i];

    file << EDGE_SE3 << " " << e.association().first << " " << e.association().second << " ";

    Vector3 t = e.data().translation();
    file << t.x() << " " << t.y() << " " << t.z() << " ";

    QuaternionReal q(e.data().linear());
    file << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " " << endl;

    OmegaPose omega = e.omega();
    file << omega.row(0)(0) << " " << omega.row(0)(1) << " " << omega.row(0)(2) << " " << omega.row(0)(3) << " " << omega.row(0)(4) << " " << omega.row(0)(5) << " " <<
            omega.row(1)(1) << " " << omega.row(1)(2) << " " << omega.row(1)(3) << " " << omega.row(1)(4) << " " << omega.row(1)(5) << " " <<
            omega.row(2)(2) << " " << omega.row(2)(3) << " " << omega.row(2)(4) << " " << omega.row(2)(5) << " " <<
            omega.row(3)(3) << " " << omega.row(3)(4) << " " << omega.row(3)(5) << " " <<
            omega.row(4)(4) << " " << omega.row(4)(5) << " " <<
            omega.row(5)(5);
    file << endl;
  }

  file.close();
  cerr << BOLDGREEN << "\t" << "File exported successfully:" << endl;
  cerr << "\t" << _vertices.size() << " Vertices SE3" << endl;
  cerr << "\t" << _edges.size() << " Edges Pose-Pose" << RESET << endl << endl;
}


void Graph::loadFromG2OFile(const std::string& filename_){
  cerr << BOLDYELLOW << "\t" << "Opening file " << filename_ << RESET << endl;
  fstream file(filename_);

  int p_idx = 0;
  int l_idx = 0;
  string line;
  while(getline(file, line)){
    stringstream ss(line);
    string element_type;

    // read datatype, if comment skip
    ss >> element_type;
    if(element_type == "#")
      continue;

    // process different elements
    if(element_type == VERTEX_3F){
/*
			int id = -1;
			PointXYZ p = PointXYZ::Zero();

			ss >> id;
			ss >> p.x() >> p.y() >> p.z();

			VertexXYZ vertex_xyz;
			vertex_xyz.setVertex(id, p, l_idx);
			addVertexXYZ(vertex_xyz);
			l_idx++;
/**/
      continue;
    }

    if(element_type == VERTEX_SE3){
      int id = -1;
      Vector3 t;
      QuaternionReal q;

      ss >> id;
      ss >> t.x() >> t.y() >> t.z();
      ss >> q.x() >> q.y() >> q.z() >> q.w();


      Pose T = Pose::Identity();
      T.linear() = q.matrix();
      T.translation() = t;

      Vertex vertex_se3;
      vertex_se3.setVertex(id, p_idx, T);
      addVertex(vertex_se3);
      p_idx++;
    }

    if(element_type == EDGE_POSE_POINT){
/*
			std::pair<int, int> IDs(-1,-1);
			int sens_id = -1;
			ss >> IDs.first >> IDs.second;
			ss >> sens_id;

			PointMeas z_land = PointMeas::Zero();
			ss >> z_land.x() >> z_land.y() >> z_land.z();

			OmegaPoint omega_land = OmegaPoint::Identity();
			ss >> omega_land.row(0)(0) >> omega_land.row(0)(1) >> omega_land.row(0)(2) >>
					omega_land.row(1)(1) >> omega_land.row(1)(2) >>
					omega_land.row(2)(2);
			for(int i = 0; i < omega_land.rows(); ++i)
				for(int j = i + 1; j < omega_land.cols(); ++j)
					omega_land(j,i) = omega_land(i,j);

			EdgePosePoint edge_pose_point;
			edge_pose_point.setEdge(IDs, sens_id, z_land, omega_land);
			addEdgePosePoint(edge_pose_point);
/**/
      continue;
    }

    if(element_type == EDGE_SE3){
      std::pair<int, int> IDs(-1,-1);
      int sens_id = -1;
      ss >> IDs.first >> IDs.second;

      Vector3 t;
      QuaternionReal q;
      ss >> t.x() >> t.y() >> t.z();
      ss >> q.x() >> q.y() >> q.z() >> q.w();

      PoseMeas odom_meas = PoseMeas::Identity();
      odom_meas.linear() = q.matrix();
      odom_meas.translation() = t;

      OmegaPose omega_odom = OmegaPose::Identity();
      ss >> 	omega_odom.row(0)(0) >> omega_odom.row(0)(1) >> omega_odom.row(0)(2) >> omega_odom.row(0)(3) >> omega_odom.row(0)(4) >> omega_odom.row(0)(5) >>
          omega_odom.row(1)(1) >> omega_odom.row(1)(2) >> omega_odom.row(1)(3) >> omega_odom.row(1)(4) >> omega_odom.row(1)(5) >>
          omega_odom.row(2)(2) >> omega_odom.row(2)(3) >> omega_odom.row(2)(4) >> omega_odom.row(2)(5) >>
          omega_odom.row(3)(3) >> omega_odom.row(3)(4) >> omega_odom.row(3)(5) >>
          omega_odom.row(4)(4) >> omega_odom.row(4)(5) >>
          omega_odom.row(5)(5);
      for(int i = 0; i < omega_odom.rows(); ++i)
        for(int j = i + 1; j < omega_odom.cols(); ++j)
          omega_odom(j,i) = omega_odom(i,j);
      Edge edge_odom;
      edge_odom.setEdge(IDs, odom_meas, omega_odom);
      addEdge(edge_odom);
    }
  }
  cerr << BOLDGREEN << "\t" << "File loaded successfully:" << endl;
  cerr << "\t" << _vertices.size() << " Vertices SE3" << endl;
  cerr << "\t" << _edges.size() << " Edges Pose-Pose" << RESET << endl << endl;
}

} /* namespace sparse */
