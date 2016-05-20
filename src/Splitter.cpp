/*
 * Splitter.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: zhonghua
 */

#include "Splitter.h"

#include <cfloat>
#include <unordered_set>
#include <algorithm>
using namespace std;

namespace masc {

Splitter::Splitter() {
  this->m_max_edge_length = FLT_MIN;
  this->m_min_edge_length = FLT_MAX;
}

void Splitter::measure(model* m) {
  for (auto i = 0; i < m->e_size; i++) {
	const auto& edge = m->edges[i];
	const auto fid1 = edge.fid[0];
	const auto fid2 = edge.fid[1];

	const auto& v1 = m->vertices[edge.vid[0]];
	const auto& v2 = m->vertices[edge.vid[1]];

	const auto edge_len = (float) ((v2.p - v1.p).norm());

	this->m_max_edge_length = std::max(this->m_max_edge_length, edge_len);
	this->m_min_edge_length = std::min(this->m_min_edge_length, edge_len);
  }

  cout << string(40, '-') << endl;
  cout << "measure result:" << endl;
  cout << "edge length = " << this->m_min_edge_length << " / "
	  << this->m_max_edge_length << endl;
}

Splitter* Splitter::createSplitter(CutHeuristic heuristic) {
  switch (heuristic) {
  case CutHeuristic::MINIMUM_PERIMETER:
	return new MinimumPerimeterSpliiter();
  case CutHeuristic::FLAT_TREE:
	return new FlatTreeSpliiter();
  case CutHeuristic::STEEPEST_EDGE:
	return new SteepestEdgeSplitter();
  case CutHeuristic::RANDOM:
	return new RandomSplitter();
  case CutHeuristic::BRUTE_FORCE:
	return new BruteForceSplitter();
  case CutHeuristic::MY_UNFOLDER_01:
	  return new MySplitter01();
  case CutHeuristic::MY_UNFOLDER_02:
	  return new MySplitter02();
  default:
	assert(false);
	break;
  }

  return nullptr;
}

//////////////////////////////////////////////////////////////////////

vector<float> MinimumPerimeterSpliiter::assignWeights(model *m,
	const Config& config) {
  vector<float> weights(m->e_size);

  for (auto i = 0; i < m->e_size; i++) {
	const auto& edge = m->edges[i];
	const auto fid1 = edge.fid[0];
	const auto fid2 = edge.fid[1];

	const auto& v1 = m->vertices[edge.vid[0]];
	const auto& v2 = m->vertices[edge.vid[1]];

	const auto edge_len = (v2.p - v1.p).norm();

	float weight = 1.0 - (edge_len - this->m_min_edge_length)
		/ (this->m_max_edge_length - this->m_min_edge_length);

	if (config.less_cuts) {
	  //cout<<"e.fa = 0"<<endl;
	  if (edge.type == 'd')
		weight = 1e10;
	  else if (edge.type == 'p')
		weight = 1e5;
	}

	weights[i] = weight;
  }

  return weights;
}

//////////////////////////////////////////////////////////////////////

vector<float> FlatTreeSpliiter::assignWeights(model *m, const Config& config) {
  vector<float> weights(m->e_size);

  // generate a random reference vector or use the given vector
  const auto c =
	  (config.use_user_vector) ?
		  config.user_vector :
		  Vector3d(mathtool::drand48(), mathtool::drand48(),
			  mathtool::drand48()).normalize();

  for (auto i = 0; i < m->e_size; i++) {
	const auto& edge = m->edges[i];
	const auto fid1 = edge.fid[0];
	const auto fid2 = edge.fid[1];

	const auto& v1 = m->vertices[edge.vid[0]];
	const auto& v2 = m->vertices[edge.vid[1]];

	const auto edge_len = (v2.p - v1.p).norm();

	float weight = fabs(c * (v2.p - v1.p)) / edge_len;

	if (config.less_cuts) {
	  if (edge.type == 'd')
		weight = 1e10;
	  else if (edge.type == 'p')
		weight = 1e5;
	}

	weights[i] = weight;
  }

  return weights;
}

//////////////////////////////////////////////////////////////////////
vector<float> SteepestEdgeSplitter::assignWeights(model* m,
	const Config& config) {
  vector<float> weights(m->e_size);

  // generate a random reference vector or use the given vector
  auto c =
	  (config.use_user_vector) ?
		  config.user_vector :
		  Vector3d(mathtool::drand48(), mathtool::drand48(),
			  mathtool::drand48()).normalize();
  if (config.use_user_vector == false) {
	for (short i = 0; i < 3; i++)
	  if (mathtool::drand48())
		c[i] = -c[i];
  }

  int top_vertex_id = -1;

  {
	// find the top vertex w.r.t to random vector c
	float max_prod = -FLT_MAX;
	for (auto i = 0; i < m->v_size; ++i) {
	  auto prod = c * Vector3d(m->vertices[i].p.get());
	  if (prod > max_prod) {
		max_prod = prod;
		top_vertex_id = i;
	  }
	}
  }

  unordered_set<int> selected_edges;

  for (auto i = 0; i < m->v_size; ++i) {
	if (i == top_vertex_id)
	  continue;

	const auto& v = m->vertices[i];

	float max_prod = -FLT_MAX;
	int steepest_eid = INT_MAX;
	int steepest_wid = INT_MAX;

	// find steepest edge
	for (auto eid : v.m_e) {
	  const auto& e = m->edges[eid];

	  if (config.less_cuts && e.type == 'd')
		continue;

	  const auto& wid = e.vid[0] == i ? e.vid[1] : e.vid[0];

	  const auto& w = m->vertices[wid];

	  const Vector3d vw = w.p - v.p;

	  const float prod = c * vw / vw.norm();

	  if (prod > max_prod) {
		max_prod = prod;
		steepest_eid = eid;
		steepest_wid = wid;
	  }
	}

	if (steepest_eid == INT_MAX) {
	  cerr << "steepest edge not found for vertex " << i << endl;
	} else {
	  selected_edges.insert(steepest_eid);
	}
  }

// assign weight
  for (auto i = 0; i < m->e_size; ++i) {
	const auto& e = m->edges[i];

	float weight = 0.0;

	if (selected_edges.count(i)) {
	  weight = 1.0;
	}

	const auto fid1 = e.fid[0];
	const auto fid2 = e.fid[1];

	//g[fid1][fid2] = g[fid2][fid1] = weight;

	weights[i] = weight;
  }

  return weights;
}

/////////////////////////////////////////////////////////////////////
vector<float> RandomSplitter::assignWeights(model *m, const Config& config) {
  vector<float> weights(m->e_size);

  for (auto i = 0; i < m->e_size; i++) {
	const auto& edge = m->edges[i];
	const auto fid1 = edge.fid[0];
	const auto fid2 = edge.fid[1];

	// random weight
	float weight = mathtool::drand48();

	weights[i] = weight;
  }

  return weights;
}

void BruteForceSplitter::init(int edges) {
  cout << "BruteForceSplitter::init edges = " << edges << endl;
  this->m_weights = vector<float>(edges);

  for (int i = 0; i < edges; ++i) {
	this->m_weights[i] = i;
  }
}

vector<float> BruteForceSplitter::assignWeights(model *m,
	const Config& config) {
  if (!this->m_inited) {
	this->init(m->e_size);
	this->m_inited = true;
  } else {
	if (!std::next_permutation(this->m_weights.begin(),
		this->m_weights.end())) {
	  cerr << "All possible permutation tried!" << endl;
	  assert(false);
	}
  }

  return this->m_weights;
}

/////////////////////////////////////////////////////////////////////
// First splitter : Greatest increase cut tree
vector<float> MySplitter01::assignWeights(model *m, const Config& config) 
{
	vector<float> weights(m->e_size); 

	// generate a random reference vector or use the given vector
	auto c =
		(config.use_user_vector) ?
		config.user_vector :
		Vector3d(mathtool::drand48(), mathtool::drand48(),
		mathtool::drand48()).normalize();
	if (config.use_user_vector == false) {
		for (short i = 0; i < 3; i++)
			if (mathtool::drand48())
				c[i] = -c[i];
	}

	int top_vertex_id = -1;

	{
		// find the top vertex w.r.t to random vector c
		float max_prod = -FLT_MAX;
		for (auto i = 0; i < m->v_size; ++i) {
			auto prod = c * Vector3d(m->vertices[i].p.get());
			if (prod > max_prod) {
				max_prod = prod;
				top_vertex_id = i;
			}
		}

	}

	unordered_set<int> selected_edges;

	for (auto i = 0; i < m->v_size; ++i) {
		if (i == top_vertex_id)
			continue;

		const auto& v = m->vertices[i];

		float max_prod = -FLT_MAX;
		int greatest_eid = INT_MAX;
		int greatest_wid = INT_MAX;

		// find greatest increase edge
		for (auto eid : v.m_e) {
			const auto& e = m->edges[eid];

			const auto& wid = e.vid[0] == i ? e.vid[1] : e.vid[0];

			const auto& w = m->vertices[wid];

			const Vector3d vw = w.p - v.p;

			const float prod = c * vw ;

			if (prod > max_prod) {
				max_prod = prod;
				greatest_eid = eid;
				greatest_wid = wid;
			}
		}

		if (greatest_eid == INT_MAX) {
			cerr << "greatest increase edge not found for vertex " << i << endl;
		}
		else {
			selected_edges.insert(greatest_eid);
		}
	}

	// assign weight
	for (auto i = 0; i < m->e_size; ++i) {
		const auto& e = m->edges[i];

		float weight = 0.0;

		if (selected_edges.count(i)) {
			weight = 1.0;
		}

		weights[i] = weight;
	}
	return weights;
}


/////////////////////////////////////////////////////////////////////
// Second splitter : Cutting the rightmost ascending edge
vector<float> MySplitter02::assignWeights(model *m, const Config& config)
{
	vector<float> weights(m->e_size);

	// generate a random reference vector or use the given vector
	auto c =
		(config.use_user_vector) ?
		config.user_vector :
		Vector3d(mathtool::drand48(), mathtool::drand48(),
		mathtool::drand48()).normalize();
	if (config.use_user_vector == false) {
		for (short i = 0; i < 3; i++)
			if (mathtool::drand48())
				c[i] = -c[i];
	}

	//initialize v_: minimal vertex of P w.r.t. c
	int bottom_vertex_id = -1;

	{
		// find the bottom vertex w.r.t to random vector c
		float min_prod = FLT_MAX;
		for (auto i = 0; i < m->v_size; ++i) {
			auto prod = c * Vector3d(m->vertices[i].p.get());
			if (prod < min_prod) {
				min_prod = prod;
				bottom_vertex_id = i;
			}
		}
	}

	//initialize b: the barycenter of P
	auto b = Vector3d();
	
	for (auto i = 0; i < m->v_size; i++) {
		
		auto vi = Vector3d(m->vertices[i].p.get());
		b += vi;
	}
	
	b = b / m->v_size;
	
	unordered_set<int> selected_edges;
	
	//find rightmost ascending edge 
	for (auto i = 0; i < m->v_size; i++) {
		
		if (i == bottom_vertex_id)
			continue;
		
		const auto& v = m->vertices[i];
		float max_det = -FLT_MAX;
		int rightmost_eid = INT_MAX;
		int rightmost_wid = INT_MAX;

		//find rightmost ascending edge incident to v
		for (auto eid : v.m_e) {
			const auto& e = m->edges[eid];

			const auto& wid = e.vid[0] == i ? e.vid[1] : e.vid[0];
			const auto& w = m->vertices[wid];

			const Vector3d vw = w.p - v.p;
			const Vector3d bv = v.p - b;

			const float prod = c * vw / vw.norm();

			const Vector3d f = bv / bv.norm();
			const Vector3d i = vw / vw.norm();


			//calculate determinant
			//det((a,b,c), (d,e,f), (g,h,i)) = aei + bfg + cdh - ceg - bdi - afh

			const float aei = c[0] * f[1] * i[2];
			const float bfg = c[1] * f[2] * i[0];
			const float cdh = c[2] * f[0] * i[0];
			const float ceg = c[2] * f[1] * i[0];
			const float bdi = c[1] * f[0] * i[2];
			const float afh = c[0] * f[2] * i[1];

			const float det = aei + bfg + cdh - ceg - bdi - afh;

			if (prod > 0 && det > max_det) {
				max_det = det;
				rightmost_eid = eid;
				rightmost_wid = wid;
			}
		}

		if (rightmost_eid == INT_MAX) {
			cerr << "steepest edge not found for vertex " << i << endl;
		}
		
		else {
			selected_edges.insert(rightmost_eid);
		}
	}
	
	// assign weight
	for (auto i = 0; i < m->e_size; ++i) {
		const auto& e = m->edges[i];
		
		float weight = 0.0;
		
		if (selected_edges.count(i)) {
			weight = 1.0;
		}
		
		weights[i] = weight;
	}
	return weights;
}


} /* namespace masc */
