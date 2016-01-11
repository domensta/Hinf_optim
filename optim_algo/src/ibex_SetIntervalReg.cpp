//============================================================================
//                                  I B E X                                   
// File        : ibex_Set.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 13 juil. 2014
//============================================================================

#include <stack>
#include <fstream>
#include "ibex_SetIntervalReg.h"

using namespace std;

namespace ibex {



SetIntervalReg::SetIntervalReg(const IntervalVector& bounding_box, double eps, bool inner) : root(new SetNodeReg(inner? __IBEX_IN__: __IBEX_UNK__)), eps(eps), bounding_box(bounding_box) {
	

}


/*SetIntervalReg::SetIntervalReg(const char* filename) : root(NULL), eps(-1), bounding_box(1) {
	load(filename);
}*/

SetIntervalReg::SetIntervalReg(const SetIntervalReg& set) : eps(set.eps), bounding_box(set.bounding_box)
{
	root = set.root->copy(NULL);
}

SetIntervalReg::SetIntervalReg(const char *filename): root(NULL), eps(-1), bounding_box(1) {
	load(filename);

}


bool SetIntervalReg::is_empty() const {
	return root==NULL;
}

bool SetIntervalReg::is_reg() const {
	return true;
}

void SetIntervalReg::contract(Sep& sep) {
	root->changeStat(__IBEX_IN__,__IBEX_UNK__);
	root->cleave(bounding_box, sep,eps);
    root->gather();
}

void SetIntervalReg::contract(Ctc& Ctcin, Ctc& Ctcout) {
	root->changeStat(__IBEX_IN__,__IBEX_UNK__);
	root->cleave(bounding_box, Ctcin, Ctcout,eps);
    root->gather();
}

void SetIntervalReg::contract(IntervalVector * box,int type) {
	vector<SetNodeReg*> neigh;
	findNeighbor(*box,&neigh);
	vector<IntervalVector> neighbox;
	double min,max;
	while(!neigh.empty()) {
		if(neigh.back()->status == type) {neighbox.push_back(findNodeBox(neigh.back()));}
		neigh.pop_back();
	}
	if(neighbox.empty()) {
		*box = IntervalVector::empty(box->size());
		return;
	}
	for(int i = 0;i<box->size();i++) {
		double min((*box)[i].ub()),max((*box)[i].lb());
		for(int j =0;j<neighbox.size();j++) {
			min = (min<neighbox.at(j)[i].lb()? min : neighbox.at(j)[i].lb());
			max = (max>neighbox.at(j)[i].ub()? max : neighbox.at(j)[i].ub());
		}
		if(min>(*box)[i].lb()) {(*box)[i] = Interval(min,(*box)[i].ub());}
		if(max<(*box)[i].ub()) {(*box)[i] = Interval((*box)[i].lb(),max);}
	}
}

void SetIntervalReg::operator&=(const SetIntervalReg& set) {
	if(set.is_reg()) {
		const SetIntervalReg& setreg = (const SetIntervalReg&) set;
		assert(bounding_box == setreg.bounding_box); // root must be the same as set interval are regular
		root->oper(setreg.root,true);
		root->gather();
	}
}


void SetIntervalReg::interBox(const IntervalVector& box,NodeType valin,NodeType valout) {
	root->operator_ir(bounding_box,box,valin,valout,true,eps);
	root->gather();
}

void SetIntervalReg::unionBox(const IntervalVector& box,NodeType valin,NodeType valout) {
	root->operator_ir(bounding_box,box,valin,valout,false,eps);
    root->gather();
}


void SetIntervalReg::operator|=(const SetIntervalReg& set) {
	if(set.is_reg()) {
		const SetIntervalReg& setreg = (const SetIntervalReg&) set;
		assert(bounding_box == setreg.bounding_box); // root must be the same as set interval are regular
		root->oper(setreg.root,false);
		root->gather();
	}
}

void SetIntervalReg::findNeighbor(const IntervalVector& box, vector<SetNodeReg*> * neigh) {
	root->findNeighbor(bounding_box,box,neigh);
}

void SetIntervalReg::findNeighbor(const IntervalVector& box, vector<SetNodeReg*> * neigh, vector<IntervalVector> * neighbox) {
	root->findNeighbor(bounding_box,box,neigh,neighbox);
}

void SetIntervalReg::findNeighbor(SetNodeReg* node, vector<SetNodeReg*> * neigh) {
    IntervalVector box(findNodeBox(node));
	root->findNeighbor(bounding_box,box,neigh);
}

void SetIntervalReg::findNeighbor(SetNodeReg* node, vector<SetNodeReg*> * neigh, vector<IntervalVector> * neighbox) {
    IntervalVector box(findNodeBox(node));
	root->findNeighbor(bounding_box,box,neigh,neighbox);
}

IntervalVector SetIntervalReg::findNodeBox(SetNodeReg * node) {
	IntervalVector bbox(bounding_box);
	IntervalVector box(bounding_box);
	node->findBox(bbox,&box);
	return box;
}

bool SetIntervalReg::contains(const IntervalVector& box,const int& type) const {
	bool flag(true);
	root->compare(bounding_box,box,type,&flag,0);
	return flag;
}

bool SetIntervalReg::intersects(const IntervalVector& box,const int& type) const{
	bool flag(false);
	root->compare(bounding_box,box,type,&flag,1);
	return flag;
}

bool SetIntervalReg::is_disjoint(const IntervalVector& box,const int& type) const{
	bool flag(true);
	root->compare(bounding_box,box,type,&flag,2);
	return flag;
}

vector<SetNodeReg*> SetIntervalReg::getLeaf() const{
	vector<SetNodeReg*> vecLeaf;
	root->getLeaf(&vecLeaf);
	return vecLeaf;
}

void SetIntervalReg::getLeaf(vector<SetNodeReg*> * node_vect,vector<IntervalVector> * node_box) {
	root->getLeaf(node_vect,node_box,bounding_box);
}

void SetIntervalReg::gather() {
	root->gather();
}

void SetIntervalReg::gather(int type) {
	root->gather(type);
}

void SetIntervalReg::visit_leaves(SetNodeReg::leaf_func func) const {
	root->visit_leaves(func, bounding_box);
}

/*std::ostream& operator<<(std::ostream& os, const SetIntervalReg& set) {

	set.root->print(os,set.bounding_box, 0);
	return os;
}*/

void SetIntervalReg::save(const char* filename) {
	std::stack<SetNodeReg*> s;

	s.push(root);

	fstream os;
	os.open(filename, ios::out | ios::trunc | ios::binary);

	os.write((char*) &eps, sizeof(double));

	int n=bounding_box.size();
	os.write((char*) &n, sizeof(int));

	for (int i=0; i<bounding_box.size(); i++) {
		double d; // to store double values
		d=bounding_box[i].lb();
		os.write((char*) &d,sizeof(double));
		d=bounding_box[i].ub();
		os.write((char*) &d,sizeof(double));
	}

	while (!s.empty()) {
		SetNodeReg* node=s.top();
		s.pop();
		if (node->is_leaf()) {
			int no_var=-1; // to store "-1" (means: leaf)
			os.write((char*) &no_var, sizeof(int));
			os.write((char*) &node->status, sizeof(NodeType));
		}
		else {
			os.write((char*) &node->var, sizeof(int));
			os.write((char*) &node->pt, sizeof(double));
			s.push(node->right);
			s.push(node->left);
		}
	}
	os.close();
}

void SetIntervalReg::load(const char* filename) {

	std::ifstream is;
	is.open(filename, ios::in | ios::binary);

	is.read((char*) &eps, sizeof(double));
	//cout << "eps=" << eps << endl;

	unsigned int n;
	is.read((char*) &n, sizeof(int));
	//cout << "n=" << n << endl;

	bounding_box.resize(n);

	for (int i=0; i<bounding_box.size(); i++) {
		double lb,ub;
		is.read((char*) &lb, sizeof(double));
		is.read((char*) &ub, sizeof(double));
		bounding_box[i]=Interval(lb,ub);
	}
	//cout << "bounding box=" << bounding_box << endl;

	int var;
	is.read((char*) &var, sizeof(int));

	double pt;
	NodeType status;

	if (var==-1) {
		is.read((char*) &status, sizeof(NodeType));
		root = new SetNodeReg(status);
		is.close();
		return;
	}

	is.read((char*) &pt, sizeof(double));
	std::stack<SetNodeReg*> s;
	root = new SetNodeReg(__IBEX_UNK__,NULL,var, pt); // left and right are both set to NULL temporarily
	s.push(root);

	SetNodeReg *node;
	SetNodeReg* subnode;

	while (!s.empty()) {

		assert(!s.top()->is_leaf());

		node = s.top();

		// =============== backtrack ======================
		if (node->left && node->right) {
			s.pop();
			continue;
		}

		is.read((char*) &var, sizeof(int));

		if (var==-1) {
			is.read((char*) &status, sizeof(NodeType));
			subnode = new SetNodeReg(status,node,(node->var+1)%n,pt);
		} else {
			is.read((char*) &pt, sizeof(double));
			subnode  =new SetNodeReg(__IBEX_UNK__,node,var,pt); // left and right are both set to NULL temporarily
		}

		if (node->left==NULL) node->left=subnode;
		else {
			assert(node->right==NULL);
			node->right=subnode;
		}

		if (var!=-1)
			s.push(subnode);
	}

	is.close();
}


/*namespace {


class NodeAndDist : public Backtrackable {
public:
	NodeAndDist() : node(NULL), dist(-1) { }

	NodeAndDist(SetNodeReg* _node) : node(_node), dist(-1) { }*/

	/**
	 * calculate the square of the distance to pt
	 * for the box of the current cell (box given in argument)
	 */
	/*void set_dist(const IntervalVector& box, const Vector pt) {
		assert(box.size()==pt.size());

		Interval d=Interval::ZERO;
		for (int i=0; i<pt.size(); i++) {
			d += sqr(box[i]-pt[i]);
		}
		dist=d.lb();
	}

	virtual std::pair<Backtrackable*,Backtrackable*> down() {
		assert(!node->is_leaf());

		SetBisect& b=*((SetBisect*) node);
		return std::pair<NodeAndDist*,NodeAndDist*>(new NodeAndDist(b.left),
													  new NodeAndDist(b.right));
	}

	SetNodeReg* node;
	double dist;
};*/

/**
 * Cell heap where the criterion is the distance to "pt"
 */
/*class CellHeapDist : public CellHeap {
public:

	virtual double cost(const Cell& c) const {
		return c.get<NodeAndDist>().dist;
	}
};
}

double SetIntervalReg::dist(const Vector& pt, bool inside) const {
	CellHeapDist heap;

	//int count=0; // for stats

	Cell* root_cell =new Cell(bounding_box);
	root_cell->add<NodeAndDist>();
	root_cell->get<NodeAndDist>().node = root;
	root_cell->get<NodeAndDist>().set_dist(bounding_box,pt);
	//count++;

	heap.push(root_cell);

	double lb = POS_INFINITY;

	while (!heap.empty()) {

		Cell* c = heap.pop();

		SetNodeReg* node = c->get<NodeAndDist>().node;

		assert(node!=NULL);

		if (node->status==(inside? __IBEX_IN__ : __IBEX_OUT__)) {
			double d=c->get<NodeAndDist>().dist;
			if (d<lb) {
				lb=d;
				heap.contract_heap(lb);
			}
		} else if (!node->is_leaf() && (    (inside && possibly_contains_in(node->status))
		                                || (!inside && possibly_contains_out(node->status)))) {
			SetBisect& b= *((SetBisect*) node);

			IntervalVector left=b.left_box(c->box);
			IntervalVector right=b.right_box(c->box);

			std::pair<Cell*,Cell*> p=c->bisect(left,right);

			p.first->get<NodeAndDist>().set_dist(left,pt);
			//count++;
			if (p.first->get<NodeAndDist>().dist<=lb) heap.push(p.first);

			p.second->get<NodeAndDist>().set_dist(right,pt);
			//count++;
			if (p.second->get<NodeAndDist>().dist<=lb) heap.push(p.second);
		}
		delete c;
	}
	//cout << " number of times a distance to a box has been computed: " << count << endl;
	return ::sqrt(lb);
}*/



} // namespace ibex
