/*
 * Graph.h
 *
 *  Created on: 06/mar/2017
 *      Author: istin
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "defs.h"
#include "Vertex.hpp"
#include "Edge.hpp"

namespace optimizer {

	typedef Vertex<Pose> VertexSE3;
	typedef Vertex<PointXYZ> VertexXYZ;
	typedef Edge<PointMeas, OmegaPoint> EdgePosePoint;
	typedef Edge<PoseMeas, OmegaPose> EdgePosePose;


	class Graph {
	public:
		Graph();
		Graph(const Graph& graph);
		Graph(const std::string& path_to_graph_);
		virtual ~Graph();

		void loadFromG2OFile(const std::string& filename_);
		void addVertexSE3(const VertexSE3& vertex_);
		void addVertexXYZ(const VertexXYZ& vertex_);
		void addEdgePosePoint(const EdgePosePoint& edge_);
		void addEdgeOdom(const EdgePosePose& edge_);
		void updateVerticesSE3(const std::vector<VertexSE3>& new_SE3_vertices_);
		void updateVerticesXYZ(const std::vector<VertexXYZ>& new_XYZ_vertices_);

		inline const std::vector<VertexSE3>& verticesSE3(void) const {return _vertices_SE3;};
		inline const std::vector<VertexXYZ>& verticesXYZ(void) const {return _vertices_XYZ;};
		inline const std::vector<EdgePosePose>& edgesPosePose(void) const {return _edges_pose_pose;};
		inline const std::vector<EdgePosePoint>& edgesPosePoint(void) const {return _edges_pose_point;};

		inline const int numSE3Vertices(void) const {return _vertices_SE3.size();};
		inline const int numXYZVertices(void) const {return _vertices_XYZ.size();};
		inline const int graphSize(void) const {return _vertices_SE3.size() * _vertices_XYZ.size();};
		inline const int numPosePoseEdges(void) const {return _edges_pose_pose.size();};
		inline const int numPosePointEdges(void) const {return _edges_pose_point.size();};

	private:
		std::vector<VertexSE3> _vertices_SE3;
		std::vector<VertexXYZ> _vertices_XYZ;

		std::vector<EdgePosePoint> _edges_pose_point;
		std::vector<EdgePosePose> _edges_pose_pose;

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	};
}
#endif /* GRAPH_H_ */
