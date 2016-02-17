//============================================================================
//                                  I B E X                                   
// File        : ibex_SetNode.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 13 juil. 2014
//============================================================================

#include "ibex_SetNodeReg.h"


using namespace std;

// =========== shortcuts ==================
#define IN         __IBEX_IN__
#define OUT        __IBEX_OUT__
#define UNK        __IBEX_UNK__
#define UNK_IN     __IBEX_UNK_IN__
#define UNK_OUT    __IBEX_UNK_OUT__
#define UNK_IN_OUT __IBEX_UNK_IN_OUT__
#define IN_TMP     __IBEX_IN_TMP__
// ========================================

namespace ibex {


char to_string(const NodeType& status) {
	switch(status) {
	case IN : return 'Y'; break;
	case OUT : return 'N'; break;
	default : return '?';
	}
}

SetNodeReg::SetNodeReg(NodeType status) : status(status), father(NULL), right(NULL), left(NULL), var(0), cflag(false) {

}

SetNodeReg::SetNodeReg(NodeType status, SetNodeReg *father): status(status), father(father), right(NULL), left(NULL), var(0), cflag(false) {

}
SetNodeReg::SetNodeReg(NodeType status,  unsigned int var): status(status), father(NULL), right(NULL), left(NULL), var(var), cflag(false) {

}

SetNodeReg::SetNodeReg(NodeType status, SetNodeReg *father,  unsigned int var): status(status), father(father), right(NULL), left(NULL), var(var), cflag(false) {

}

SetNodeReg::SetNodeReg(NodeType status, SetNodeReg *father,  unsigned int var, double pt): status(status), father(father), right(NULL), left(NULL), var(var), cflag(false), pt(pt) {

}

SetNodeReg::~SetNodeReg() {
	delete right;
	delete left;


}

bool SetNodeReg::is_leaf() const {
	if(right == NULL && left == NULL)
		return true;
	else
		return false;
}

SetNodeReg * SetNodeReg::copy( SetNodeReg * nodeFather) const {
	if(is_leaf())
		return new SetNodeReg(status,nodeFather,var,pt);
	else {
		SetNodeReg * node = new SetNodeReg(status,nodeFather,var,pt);
		node->left = left->copy(node);
		node->right = right->copy(node);
		return node;
	}
}

SetNodeReg * SetNodeReg::getRight() {
	return right;
}

SetNodeReg * SetNodeReg::getLeft() {
	return left;
}

void SetNodeReg::cut(const IntervalVector& box) {
	assert(is_leaf());
	unsigned int cutvar;
	if(father == NULL)
		cutvar = 0;
	else
		var = (father->var+1)%box.size();
	pt = box[var].mid();
	left = new SetNodeReg(status,this,(var+1)%box.size());
	right = new SetNodeReg(status,this,(var+1)%box.size());

}

void SetNodeReg::cut(const IntervalVector& box,unsigned int cutvar) {
	assert(is_leaf());
	var = cutvar;
	pt = box[var].mid();
	left = new SetNodeReg(status,this,0);
	right = new SetNodeReg(status,this,0);
}

void SetNodeReg::operator_ir(const IntervalVector& box,const IntervalVector& nodebox, NodeType val, bool op, double eps) {
	if (op)
		operator_ir(box, nodebox, val, IN,op, eps);
	else 
		operator_ir(box, nodebox, val, OUT,op, eps);
}

void SetNodeReg::operator_ir(const IntervalVector& box,const IntervalVector& nodebox, NodeType valin, NodeType valout, bool op, double eps) {
	int r = Interset(box,nodebox);
	switch (r) {
		case 2:
			if(op && valin != IN) 
				inter(valin);
			else if(!op && valin!= OUT)
				_union(valin);
			break;
		case 1:
			if(is_leaf()) {
				if( box[var].diam()>eps) {
					cut(box);
					left->operator_ir(left_box(box),nodebox,valin,valout,op,eps);
					right->operator_ir(right_box(box),nodebox,valin,valout,op,eps);
				}
				else
					status = inte_in(status,UNK);
			}
			else {
				left->operator_ir(left_box(box),nodebox,valin,valout,op,eps);
				right->operator_ir(right_box(box),nodebox,valin,valout,op,eps);
			}
			break;
		case 0:
			if(op && valout!= IN) 
				inter(valout);
			else if(!op && valout!= OUT)
				_union(valout);
			break;
	}

}

void SetNodeReg::oper(SetNodeReg* node,bool op) {
	if(is_leaf())
	{
		if(node->is_leaf()) // if node is leaf, apply its value to leaves of this
		{
			if(op)
				status = inte(status,node->status);
			else
				status = uni(status,node->status);
		}
		else {
			pt = node->pt;
			left = node->left->copy(this);
			right = node->right->copy(this);
			if(op && status != IN) 
			{
				left->inter(status);
				right->inter(status);
			}
			else if(!op && status!= OUT)
			{
				left->_union(status);
				right->_union(status);
			}
			
		}		
	}
	else {
		if(node->is_leaf()) // apply status to leaves of the subtree which have this as root
		{
			if(op && node->status != IN) 
			{
				left->inter(node->status);
				right->inter(node->status);
			}
			else if(!op && node->status!= OUT)
			{
				left->_union(node->status);
				right->_union(node->status);
			}
		}
		else
		{
			left->oper(node->left,op);
			right->oper(node->right,op);
		}
	}
}

void SetNodeReg::cleave(const IntervalVector& box, Sep& sep, const double eps) {
	if(is_leaf() && status == OUT)
		return;
	IntervalVector box1(box);
	IntervalVector box2(box);
	sep.separate(box1,box2);
    if(box1.is_empty()){
        inter(OUT);
//        cout<<"box: "<<box<<" set to OUT"<<endl;
    }

    else if(box2.is_empty()){
        _union(IN);
        cout<<"box: "<<box<<" set to IN"<<endl;
    }


	else // continu until box1 and box2 are disjoint
	{	
		if(is_leaf()) {
			if(box[var].diam()>eps) {
				cut(box);
				left->cleave(left_box(box),sep,eps);
				right->cleave(right_box(box),sep,eps);
			}
			else
				status = inte(status,UNK);
		}
		else {
			left->cleave(left_box(box),sep,eps);
			right->cleave(right_box(box),sep,eps);
		}
	}
}


void SetNodeReg::cleave(const IntervalVector& box, Ctc& Ctcin, Ctc& Ctcout, const double eps) {
	if(is_leaf() && status == OUT)
		return;
	IntervalVector box1(box);
	IntervalVector box2(box);
    //try{Ctcin.contract(box1);}
    //catch(EmptyBoxException&){inter(OUT);return;}
    Ctcin.contract(box1);
    if(box1.is_empty()) {
        inter(OUT);
        return;
    }
    //try{Ctcout.contract(box2);}
    //catch(EmptyBoxException&){_union(IN);return;}
    Ctcout.contract(box2);
    if(box2.is_empty()) {
        _union(IN);return;
    }
			
	if(is_leaf()) {
		if(box[var].diam()>eps) {
			cut(box);
			left->cleave(left_box(box),Ctcin,Ctcout,eps);
			right->cleave(right_box(box),Ctcin,Ctcout,eps);
		}
		else
			status = inte(status,UNK);
	}
	else {
		left->cleave(left_box(box),Ctcin,Ctcout,eps);
		right->cleave(right_box(box),Ctcin,Ctcout,eps);
	}
}


void SetNodeReg::inter(NodeType x_status) {
	if(is_leaf())
		status = inte(status,x_status);
	else {
		left->inter(x_status);
		right->inter(x_status);
	}

}
void SetNodeReg::_union(NodeType x_status) {
	if(is_leaf())
		status = uni(status,x_status);
	else {
		left->_union(x_status);
		right->_union(x_status);
	}
}

void SetNodeReg::changeStat(const NodeType& statp, const NodeType& statn) {
	if(is_leaf()) {
		if(status == statp)
			status = statn;
	}
	else {
		left->changeStat(statp,statn);
		right->changeStat(statp,statn);
	}
}

void SetNodeReg::gather() {
	if(is_leaf()) // this case should not append, exept if the root is a leaf
		return;
	else {
		left->gather();
		right->gather();
		if(right->is_leaf()&&left->is_leaf() && right->status == left->status) {
				status = left->status;
				delete left;
				delete right;
				left = NULL;
				right = NULL;
				return;
		}
	}
	
}

void SetNodeReg::gather(int type) {
	if(is_leaf()) // this case should not append, exept if the root is a leaf
		return;
	else {
        left->gather(type);
        right->gather(type);
        if((right->is_leaf()) && (left->is_leaf()) && (right->status == left->status) && (left->status == type) ) {
				status = left->status;
				delete left;
				delete right;
				left = NULL;
				right = NULL;
				return;
		}
	}

}

void SetNodeReg::gatherFrom() {
	gather();
	if(is_leaf() && father != NULL) {
			father->gatherFrom();
	}	
}

void SetNodeReg::gatherFrom(int type) {
	gather(type);
	if(is_leaf() && father != NULL) {
			father->gatherFrom(type);
	}	
}

void SetNodeReg::findNeighbor(const IntervalVector& box,const IntervalVector& nbox,vector<SetNodeReg*> * neigh) {
	if(is_leaf()){
		if(box.intersects(nbox)  && box!=nbox ) 
			neigh->push_back(this);
	}
    else if(box.intersects(nbox)) {
            right->findNeighbor(right_box(box),nbox,neigh);
            left->findNeighbor(left_box(box),nbox,neigh);
        }
}

void SetNodeReg::findNeighbor(const IntervalVector& box,const IntervalVector& nbox,vector<SetNodeReg*> * neigh, vector<IntervalVector> * neigbox) {
	if(is_leaf()){
		if(box.intersects(nbox) && box!=nbox ) {
			neigh->push_back(this);
			neigbox->push_back(box);}
	}
    else if(box.intersects(nbox)) {
            right->findNeighbor(right_box(box),nbox,neigh,neigbox);
            left->findNeighbor(left_box(box),nbox,neigh,neigbox);
        }
}
void SetNodeReg::findBox(const IntervalVector& bbox,IntervalVector * box) {
	if(father != NULL) {
		father->findBox(bbox,box);
		if(this == father->right) {
			*box = father->right_box(*box);
		}
		else {
			*box = father->left_box(*box);
		}
	}
}

void SetNodeReg::compare(const IntervalVector& root, const IntervalVector& box,const int& type, bool * flag,unsigned int op) {
	if(is_leaf()) { 
		//cout<<"leaf box: "<<root<<"of status: "<<status<<" box: "<<box<<endl;
		switch(op){
			case 0: 
				if(status != type){*flag = false;}
				break;
			case 1: 
				if(status == type){*flag = true;}
				break;
			case 2: 
				if(status == type){*flag = false;}
				break;
			}
	}
	else if(((op==0||op==2)&& (*flag)) || (op ==1 && !(*flag))) {
		if(root.intersects(right_box(box)))
		{right->compare(right_box(root),box,type,flag,op);}
		if(root.intersects(left_box(box)))
		{left->compare(left_box(root),box,type,flag,op);}
	}
	return;


}
void SetNodeReg::getLeaf(vector<SetNodeReg*> * vecLeaf) {
	if(is_leaf())
		vecLeaf->push_back(this);
	else {
		right->getLeaf(vecLeaf);
		left->getLeaf(vecLeaf);
	}
}

void SetNodeReg::getLeaf(vector<SetNodeReg*> * vecLeaf, vector<IntervalVector> * vecBox, const IntervalVector& curbox) {
	if(is_leaf()) {
		vecLeaf->push_back(this);
		vecBox->push_back(curbox);
	}
	else {
		right->getLeaf(vecLeaf,vecBox,right_box(curbox));
		left->getLeaf(vecLeaf,vecBox,left_box(curbox));
	}
}

void SetNodeReg::visit_leaves(leaf_func func, const IntervalVector& nodebox) const {
	if(is_leaf())
		func(nodebox, status==IN? YES : (status==OUT? NO : MAYBE));
	else {
		left->visit_leaves(func, left_box(nodebox));
		right->visit_leaves(func, right_box(nodebox));
	}
}

IntervalVector SetNodeReg::left_box(const IntervalVector& nodebox) const {
	IntervalVector leftbox(nodebox);
	leftbox[var] =Interval(nodebox[var].lb(),pt);
	return leftbox;
}

IntervalVector SetNodeReg::right_box(const IntervalVector& nodebox) const {
	IntervalVector rightbox(nodebox);
	rightbox[var]=Interval(pt,nodebox[var].ub());
	return rightbox;
}

void SetNodeReg::print(std::ostream& os, const IntervalVector& nodebox, int shift) const {
	if(is_leaf()) {
		for (unsigned int i=0; i<shift; i++) os << ' ';
		os << "* " << nodebox << endl;
	}
	else {
		for (int i=0; i<shift; i++) os << ' ';
		os  << nodebox << " " << to_string(status) << endl;
		left->print(os, left_box(nodebox), shift+2);
		right->print(os, right_box(nodebox), shift+2);
	}
}

int Interset(IntervalVector x, IntervalVector y) {
assert(x.size() == y.size());
int res = 2;
int i;
for(i = 0;i<x.size();i++) {
	if (!x[i].is_subset(y[i])) {
		res = 1;
		break;}
	}
while(i<x.size()) {
    if(!x.intersects(y)) {
		res = 0;
		return res;
	}
	i++;}
return res;}
bool is_point(IntervalVector v) {
    for(int i=0;i<v.size();i++) {
        if(!(v[i].lb()==v[i].ub())) return false;
    }
    return true;
}


} // namespace ibex
