//
//  main.cpp
//  AMP and SM1 ES
//
//  Created by Raffaele Marino on 2018-08-25.
//  Copyright © 2018 Raffaele Marino. All rights reserved.
//


#include "Header.h"
#include "Random.hpp"


/*class vertex*/
class Vertex{
public:
    /*public members*/
    Vertex(void){
        _vertex=0;
        _message_i=complex<long double>(0.,0.);
        _degree_vertex_i=0;
        _I_am_fixed=true;
    }
    
    ~Vertex(){};
    
    /*public variables*/
    friend ostream& operator<<(ostream &out, Vertex &V){
        out<<V._vertex<<" "<<real(V._message_i)<<endl;
        return out;
    }
    /*N.B.: real(y) is the message a time t; imag(y) is the message a time t+1*/
    vector<complex<long double> > _message_i_to_j; /*messages i->j*/
    complex<long double> _message_i; /*message vertex i*/
    long unsigned int _vertex; /*vertex i*/
    long unsigned int _degree_vertex_i; /*degree vertex i*/
    bool _I_am_fixed;
};

/*class graph*/
class Graph{
public:
    
    /*public members*/
    Graph(int argc, const char * argv[]){
       _N=(long unsigned int)atoi(argv[argc-3]);
        _alpha=(long double)atof(argv[argc-2]);
        _HC=(long double)atof(argv[argc-1]);
        _kQ.resize(_N); //vector useful for computing jngb
        _flag_SM0=false;
        _flag_SM1=false;
        _c_success=0;
        _m_success=0;
        cout<<_HC<<endl;
        _update_Graph();
    };
    
    Graph(long unsigned int &N, long unsigned int &M){
        _N=N;
        _M=M;
        _HC=0;
        _kQ.resize(_N); //vector useful for computing jngb
        _flag_SM0=false;
        _c_success=0;
        _m_success=0;
        _alpha=(long double)(_M*_N);
        _update_Graph();
    };
    
    ~Graph(){};
    
    /*Matula values*/
    
    const long unsigned int Matulavalue(){
       
        long double prob=(long double) _N;
        long unsigned int n=_N;
        long unsigned int matulaval=0;
        for (long unsigned int k_ = 2;k_<_N;++k_){
            --n;
            prob = prob*((long double) n)/((long double) k_);
            if(prob*pow(_alpha, (long double)((k_*(k_-1)))/2.)<0.99){
                matulaval=k_- 1;
                break;
            }
            
        };
        cout<<"The theoretical value of the maximum clique: "<<matulaval<<endl;
        
        return matulaval;
    }
    
    
    
    
    /*public member adj matrix */
    inline bool ADJ(int &i, int &j){return _adjmatrix[i][j];};
    inline bool ADJ(unsigned int &i, unsigned int &j){return _adjmatrix[i][j];};
    inline bool ADJ(long unsigned int &i, long unsigned int &j){return _adjmatrix[i][j];};
    
    /*functions that return values + or - 1*/
    inline const long double _cast_bool_to_int(long unsigned int &i, long unsigned int &j){
        return (_adjmatrix[i][j]) ? 1.:-1.;
    }
    inline const long double _cast_bool_to_int(unsigned int &i, unsigned int &j){
        return (_adjmatrix[i][j]) ? 1.:-1.;
    }
    inline const long double _cast_bool_to_int(int &i, int &j){
        return (_adjmatrix[i][j]) ? 1.:-1.;
    }
    
    /*build a E-R random graph with edge density alpha, and with a HC*/
    void build_graph(vector<Vertex *> &V){
        cout<<"This is the graph number: "<<++_counter<<endl;
        _clean_ngbs(V);
        for (long unsigned int i=0; i<_N; ++i) {
            long double rn=0.;
            for (long unsigned int j=i+1; j<_N; ++j) {
                rn=(long double)log(MyRandformessages.random_number());
                V[i]->_message_i_to_j[j]=complex<long double>(rn,rn);
                rn=(long double)log(MyRandformessages.random_number());
                V[j]->_message_i_to_j[i]=complex<long double>(rn,rn);
                if(MyRand.random_number()>_alpha){
                    _adjmatrix[i][j]=_adjmatrix[j][i]=true;
                    _ngbs[i].push_back(j);
                    _ngbs[j].push_back(i);
                    V[i]->_degree_vertex_i++;
                    V[j]->_degree_vertex_i++;
                }else{
                    _adjmatrix[i][j]=_adjmatrix[j][i]=false;
                }
            }
        }
        /*build the complete subgraph between 0 and HC*/
        for (long unsigned int i=0; i<_HC; ++i) {
            for (long unsigned int j=i+1; j<_HC; ++j) {
                if(!_adjmatrix[i][j]){
                    _adjmatrix[i][j]=_adjmatrix[j][i]=true;
                    V[i]->_degree_vertex_i++;
                    V[j]->_degree_vertex_i++;
                    for (list<long unsigned int>::iterator _it=_ngbs[i].begin(); _it!=_ngbs[i].end(); _it++) {
                        if(j<*_it){
                            _ngbs[i].insert(_it,j);
                            break;
                        }
                    }
                    for (list<long unsigned int>::iterator _it=_ngbs[j].begin(); _it!=_ngbs[j].end(); ++_it) {
                        if(i<*_it){
                            _ngbs[j].insert(_it,i);
                            break;
                        }
                    }
                }
            }
        }
        /*I check that the hidden clique is correctly planted*/
        for (unsigned int i=0; i<_HC; ++i) {
            for (unsigned int j=i+1; j<_HC; ++j) {
                if(!_adjmatrix[i][j])exit(-1);
                if(!_adjmatrix[j][i])exit(-1);
            }
        }
    };
    
    /*build the frontier for the greedy search*/
    long unsigned int frontier_site_i_given_clique(long unsigned int x, long unsigned int y){
        unsigned long i=0, frontier_site_i=0;
        for (i=0; i<_N;++i) {
            _kQ[i]=0;
        }
        i=0;
        if(!_adjmatrix[y][x]){cout<<"PRENDI MALE I NODI"<<endl; exit(-1);}
        _kQ[i] = x;
        _kQ[++i]=y;
        for(list<long unsigned int>::iterator it=_ngbs[x].begin(); it!=_ngbs[x].end(); it++){
            if(_adjmatrix[y][*it])_kQ[++i]=*it;
        }
        
        frontier_site_i=i;
        return frontier_site_i;
    }
    
    
    /*build the frontier for the greedy search*/
    long unsigned int frontier_site_i_given_clique(long unsigned int x){
        unsigned long i=0, frontier_site_i=0;
        for (i=0; i<_N;++i) {
            _kQ[i]=0;
        }
        i=0;
        _kQ[i] = x;
        for(list<long unsigned int>::iterator it=_ngbs[x].begin(); it!=_ngbs[x].end(); it++){
            _kQ[++i]=*it;
        }
        
        frontier_site_i=_ngbs[x].size();
        return frontier_site_i;
    }
    
    /*greedy algorithm*/
    vector<long unsigned int> greedyloop(long unsigned int &jngb, int s=1){
        vector<long unsigned int> cliquefound;
        int imin, maxlinks, kk, ovlp;
        for (int npruned=s;npruned < jngb;++npruned){
            imin = (int) _N;maxlinks = 0;kk = 0;
            for (int i=npruned;i<jngb;++i){  /*calculate ovlp of site i   with rest of the clique horizon*/
                ovlp = 0;
                for (int j=0;j<jngb;++j)
                    if (_adjmatrix[_kQ[i]][_kQ[j]])    ++ovlp;
                if (imin>ovlp) imin = ovlp;
                if (maxlinks<ovlp){maxlinks = ovlp;kk=i;
                }
            }
            
            long unsigned int temp = _kQ[kk];
            _kQ[kk] = _kQ[npruned];
            _kQ[npruned] = temp;
            
            int ii = npruned+1;
            for (int i=npruned+1;i<jngb;++i){
                if (_adjmatrix[_kQ[i]][_kQ[npruned]]){
                    _kQ[ii] = _kQ[i];++ii;
                }
            }
            jngb = ii;
            
            if (jngb == npruned+1) {
                /*check if it is a clique*/
                for (long unsigned int i=0;i<=npruned;++i) {cliquefound.push_back(_kQ[i]);}
                sort(cliquefound.begin(), cliquefound.end());
                for (long unsigned int i=0; i<cliquefound.size(); ++i) {
                    for (long unsigned int j=i+1; j<cliquefound.size(); ++j) {
                        if(!_adjmatrix[cliquefound[i]][cliquefound[j]]){
                            cout<<_adjmatrix[cliquefound[i]][cliquefound[j]]<<" "<<cliquefound[i]<<" "<<cliquefound[j]<<endl;
                            cout<<"clique not found :("<<endl;
                            exit(-1);/*stop completely the program if it is not a clique*/
                        }
                    }
                    
                }
            }
        }
        return cliquefound;
    }
    
    /*cleaning algorithm*/
    void cleaning(vector<long unsigned int> &cliquefound){
        vector<bool> flag_cl;
        vector<complex<long unsigned int> > vec_rec;
        vector<complex<long unsigned int> > cpvec_rec;
        vector<complex<long unsigned int> > vec_rec_cl;
        vector<long unsigned int> var;
        vector<long unsigned int> clcounter;
        vector<long unsigned int> cl_old;
        cl_old=cliquefound;
        t__=0;
        for (unsigned long t=0; t<20; ++t) {
            cout<<"Iterazione: "<<t<<" "<<cliquefound.size()<<endl;
            if(cliquefound.size()==_HC)break; // break
            t__=t;
            flag_cl.clear();
            flag_cl.resize(_N);
            for(unsigned long i=0; i<cliquefound.size(); ++i){
                flag_cl[cliquefound[i]]=true;
            }
            vec_rec.clear();
            cpvec_rec.clear();
            vec_rec_cl.clear();
            clcounter.clear();
            
            clcounter.resize(cliquefound.size());
            vec_rec.resize(_N);
            var.clear();
            var.resize(_N);

            for (unsigned long i=0; i<_N; ++i) {
                 // count the number of connections that site i has with all elements in the set C
                if(!flag_cl[i])for (unsigned long j=0; j<cliquefound.size(); ++j) {
                    if(_adjmatrix[i][cliquefound[j]]){
                        ++clcounter[j];
                        ++var[i];
                    }
                }
                vec_rec[i]=complex<unsigned long>(var[i],i);
                
            }
            // take the higest site and put them into the set C
            sort(vec_rec.begin(), vec_rec.end(), _complex_greater_pred<unsigned long>());
            unsigned long old=cliquefound.size();
            unsigned long c=1;
            cliquefound.push_back((int)imag(vec_rec[c-1]));
            while (real(vec_rec[c-1])==real(vec_rec[c])) {
                cliquefound.push_back((int)imag(vec_rec[c]));
                ++c;
            }
            vec_rec_cl.resize(cliquefound.size());
            clcounter.clear();
            clcounter.resize(cliquefound.size());
            // count the  number of links inside set C and pick the ones that have the higest degree
            for (unsigned long i=0; i<cliquefound.size(); ++i) {
                for (unsigned long j=0; j<cliquefound.size(); ++j) {
                    if(_adjmatrix[cliquefound[i]][cliquefound[j]]){
                        ++clcounter[i];
                    }
                }
                
                vec_rec_cl[i]=complex<unsigned long>(clcounter[i],cliquefound[i]);
                
            }
            
            sort(vec_rec_cl.begin(), vec_rec_cl.end(), _complex_greater_pred<unsigned long>());

            cliquefound.clear();
            c=1;
            cliquefound.push_back((int)imag(vec_rec_cl[c-1]));
            while (real(vec_rec_cl[c])>old and c<vec_rec_cl.size()) {
                cliquefound.push_back((int)imag(vec_rec_cl[c]));
                ++c;
                
            }
            
            
            
        }
        /*check if it is a clique*/
        for (unsigned long i=0; i<cliquefound.size(); ++i) {
            for (unsigned long j=i+1; j<cliquefound.size(); ++j) {
                if(!_adjmatrix[cliquefound[i]][cliquefound[j]]){
                    cout<<"ERROR!!!!!!"<<cliquefound[i]<<" "<<cliquefound[j]<<endl;
                    cliquefound.clear();
                    cliquefound=cl_old;
                    break;
                }
            }
        }
        
        sort(cliquefound.begin(), cliquefound.end());
        for (unsigned long i=0; i<cliquefound.size(); ++i) { // print the clique
            cout<<i<<" "<<cliquefound[i]<<endl;
        }
  
        if(cliquefound.size()>cl_old.size()){
        vec_rec_final.push_back(cliquefound.size());
        }else{
        vec_rec_final.push_back(cl_old.size());
        }
     
    }
    /*SM1 ES*/
    void SM_2_and_stop(vector<Vertex *> &V){
        _flag_SM1=true;
        unsigned long rn;
        for (unsigned long i=0; i<_N; ++i) {// randomise the elements in vector V
            rn=rand()%_N;
            swap(V[i], V[rn]);
        }
        long unsigned int jngb=0, iter=0;
        const long unsigned int ml=Matulavalue();//compute the value of the staircase associated to N
        vector<long unsigned int> cl;
        bool flag_sm2=false;
        for (long unsigned int i=0; i<_N; ++i) {
            for (long unsigned int j=i+1; j<_N; ++j) { // double for because we are looking for pairs
                if(_adjmatrix[V[i]->_vertex][V[j]->_vertex]){ // check is exist the edge between vertecies
                    cout<<*V[i]<<" "<<*V[j]<<" "<<_adjmatrix[V[i]->_vertex][V[j]->_vertex]<<endl;
                    jngb=frontier_site_i_given_clique(V[i]->_vertex,V[j]->_vertex);  // build the frontier
                    cl.clear();
                    cl=greedyloop(jngb, 2); // greedy algorithm
                    cout<<"size SM("<<2<<"): "<<cl.size()<<endl;
                    sort(cl.begin(), cl.end());
                    iter=i+1;
                    if(cl.size()>(ml+2)){
                        flag_sm2=true;
                        break;
                    }
                }
            }
            if(flag_sm2)break;
        }
        
        
       
        unsigned long int _counter_HC_sites=0;
        for (unsigned long i=0; i<cl.size(); ++i) {
            if(cl[i]<_HC)++_counter_HC_sites;
        }
        _frac_clique_found_SM1.push_back(_counter_HC_sites);
        _size_clique_hc_sites.push_back(complex<unsigned long>((unsigned long)_counter_HC_sites,cl.size()));
        _mi.push_back(iter);
        if(cl.size()>ml+2){
            t__=0;
            cleaning(cl); // cleaning algorithm is applied on the clique found
        }else{
            t__=20;
        }
        _iteration_cleaning.push_back(complex<unsigned long>(t__,cl.size()));
        if(cl.size()==_HC)++_c_success;
        sort(V.begin(), V.end(), _pred_sort_smaller_vertex());
    
    }
    
    /*SM1 ES*/
    void SM_1_and_stop(vector<Vertex *> &V){
        _flag_SM1=true;
        unsigned long rn;
        for (unsigned long i=0; i<_N; ++i) {// randomise the elements in vector V
            rn=rand()%_N;
            swap(V[i], V[rn]);
        }
        long unsigned int jngb=0, iter=0;
        const long unsigned int ml=Matulavalue(); //compute the value of the staircase associated to N
        vector<long unsigned int> cl;
        for (long unsigned int i=0; i<_N; ++i) {
            cout<<*V[i];
            jngb=frontier_site_i_given_clique(V[i]->_vertex); // build the frontier
            cl.clear();
            cl=greedyloop(jngb, 1); // find the clique using greey
            cout<<"size SM("<<1<<"): "<<cl.size()<<endl; // print the size of the clique
            sort(cl.begin(), cl.end());
            iter=i+1;
            if(cl.size()>(ml+2))break; // if the size is bigger than ml+2 stop
        }
        unsigned long int _counter_HC_sites=0;
        for (unsigned long i=0; i<cl.size(); ++i) {
            if(cl[i]<_HC)++_counter_HC_sites;
        }
        _frac_clique_found_SM1.push_back(_counter_HC_sites);
        _size_clique_hc_sites.push_back(complex<unsigned long>((unsigned long)_counter_HC_sites,cl.size()));
        _mi.push_back(iter);
        if(cl.size()>ml+2){
            t__=0;
            cleaning(cl); // cleaning algorithm is applied on the clique found
        }else{
            t__=20;
        }
        _iteration_cleaning.push_back(complex<unsigned long>(t__,cl.size()));
        if(cl.size()==_HC)++_c_success;
        sort(V.begin(), V.end(), _pred_sort_smaller_vertex());
    }
    
    /*iterative procedure for BP equations for SM0*/
    /*Warning function not optimizsed*/
    void _iterative_procedure_for_SM0(vector<Vertex *> &V){
        _flag_SM0=true;
        long double frac_HC=0.;
        long unsigned int max=0.;
        for (unsigned int t=0; t<t_max; ++t) {
            max=0.;
            frac_HC=0;
            cout<<"iteration: "<<t<<endl;
            if(_update_messages(V)){
                cout<<"I found a convergence"<<endl;
                cout<<"iteration: "<<t<<endl;
                sort(V.begin(), V.end(), _pred_sort_greater()); /*sort for fraction recovery*/
                for (long unsigned int i=0; i<_HC; ++i) {
                    cout<<*V[i];
                    if(V[i]->_vertex<_HC){
                        frac_HC++;
                    }
                }
                if(frac_HC!=_HC){
                    frac_HC=0;
                long unsigned int jngb=0, iter=0;
                const long unsigned int ml=Matulavalue();
                vector<long unsigned int> cl;
                for (long unsigned int i=0; i<100; ++i) {
                    // da questo punto in poi devo iniziare a scrivere il codice per implementare SM0 sui primi 100 siti
                    cout<<*V[i];
                    jngb=frontier_site_i_given_clique(V[i]->_vertex);
                    cl.clear();
                    cl=greedyloop(jngb, 1);
                    cout<<"size SM("<<0<<"): "<<cl.size()<<endl;
                    sort(cl.begin(), cl.end());
                    
                    if(V[i]->_vertex<_HC){
                        frac_HC++;
                    }
                    iter=i+1;
                    if(cl.size()>(ml+5))break;
                }
                _frac_clique_found.push_back(cl.size());
                _mi.push_back(iter);
                if(cl.size()>ml)cleaning(cl);
                if(cl.size()==_HC)++_c_success;
                cout<<"Number of sites of HC in the first 100 position is: "<<frac_HC<<endl;
                }else{
                    _frac_clique_found.push_back(_HC);
                    _m_success++;
                    _mi.push_back(0);
                    
                }
                t=t_max;
            }
        }
    }
    
    
    
    
    /*iterative procedure for BP equations*/
    void _iterative_procedure(vector<Vertex *> &V){
        long double frac_HC=0.;
        long unsigned int max=0.;
        unsigned long t_=0;
        for (unsigned int t=0; t<t_max; ++t) {
            max=0.;
            frac_HC=0;
            cout<<"iteration: "<<t<<endl;
            t_=t;
            if(_update_messages(V)){
                cout<<"I found a convergence"<<endl;
                cout<<"iteration: "<<t<<endl;
                sort(V.begin(), V.end(), _pred_sort_greater()); /*sort for fraction recovery*/

                for (long unsigned int i=0; i<_HC; ++i) {
                    cout<<*V[i];
                    if(V[i]->_vertex<_HC){
                        frac_HC++;
                    }
                }
                _frac_clique_found.push_back(frac_HC);
                if(frac_HC==_HC)_m_success++;
                t=t_max;
                sort(V.begin(), V.end(), _pred_sort_smaller_vertex());
            }
        }
        _size_clique_hc_sites_AMP.push_back(complex<unsigned long>((unsigned long)frac_HC,(unsigned long)_HC));
        _iteration_AMP.push_back(complex<unsigned long>((unsigned long)t_,(unsigned long)frac_HC));
    }
    
    
    void _ri_initialization(vector<Vertex *> &V){
        for (long unsigned int i=0; i<_N; ++i) {
            if(V[i]->_I_am_fixed)V[i]->_message_i=complex<long double>(0.,0.);
            else V[i]->_message_i=complex<long double>(-((long double)_N),-((long double)_N));
            for (long unsigned int j=i+1; j<_N; ++j) {
                V[i]->_message_i_to_j[j]=complex<long double>((long double)0.,(long double)0.);
                V[j]->_message_i_to_j[i]=V[i]->_message_i_to_j[j];
            }
        }
    }
    

    friend ostream& operator<<(ostream& out, Graph& B){ /*operator << for object graph*/
        const long double mean=B._mean(B._frac_clique_found);
        const long double mean1=B._mean(B._mi);
        const long double mean2=B._mean(B._frac_clique_found_SM1);
        const long double mean3=B._mean(B.vec_rec_final);
        if(B._flag_SM0){
        out<<B.order_Graph<<" "<<B.size_HC<<" "<<mean<<" "<<B._std(B._frac_clique_found,mean)<<" "<<mean1<<" "<<B._std(B._mi,mean1)<<" "<<mean2<<" "<<B._std(B.vec_rec_final,mean2)<<" "<<B._c_success<<" "<<B._m_success<<" "<<NGraphs<<endl;
        }else if(B._flag_SM1){
        out<<B.order_Graph<<" "<<B.size_HC<<" "<<mean<<" "<<B._std(B._frac_clique_found,mean)<<" "<<mean2<<" "<<B._std(B._frac_clique_found_SM1,mean2)<<" "<<mean1<<" "<<B._std(B._mi,mean1)<<" "<<mean3<<" "<<B._std(B.vec_rec_final,mean3)<<" "<<B._c_success<<" "<<B._m_success<<" "<<NGraphs<<endl;
        }else{
            out<<B.order_Graph<<" "<<B.size_HC<<" "<<mean<<" "<<B._std(B._frac_clique_found,mean)<<" "<<hint<<endl;
        }
        return out;
    }
    
    long unsigned int give_me_HC(){return _HC;}
    /*public variables*/
    Random MyRand; /*random number generator*/
    Random MyRandformessages; /*random number generator for messages*/
    vector<list<long unsigned int> > _ngbs; /*list of ngbs*/
     vector<complex<unsigned long> > _size_clique_hc_sites_AMP;
    vector<complex<unsigned long> > _iteration_AMP;
    vector<complex<unsigned long> > _size_clique_hc_sites;
    vector<complex<unsigned long> > _iteration_cleaning;
    vector<long double> _frac_clique_found; /*fraction of clique recovered*/
    vector<long double> _frac_clique_found_SM1; /*fraction of clique recovered*/
    vector<long double> _mi, vec_rec_final; /*number_of_sites_traveled*/
    long double _inv_sqrtN; /*inverse of SQRT N*/
    long unsigned int order_Graph; /*Graph order*/
    long unsigned int size_Graph; /*Graph size*/
    long unsigned int size_HC; /*Hidden clique size*/
    unsigned int _counter; /*counter for number of graphs*/
    long unsigned int _c_success; /*count success SM1-ES*/
    long unsigned int _m_success;/*count success AMP*/
    unsigned long t__;
    bool _flag_SM0, _flag_SM1;
    
    
    
private:
    
    /*private variables*/
    vector<complex<long double> > _vec_sort_i_mar; /**/
    vector<vector<bool> > _adjmatrix; /*Adjacency matrix*/
    vector<long unsigned int> _kQ; /* vector for greedy algorithm*/
    long double _alpha; /*edge density*/
    long double _kappa; /*constant values for messages equations*/
    long unsigned int _N; /*number of nodes*/
    long unsigned int _M; /*number of edges*/
    long unsigned int _HC; /*Hidden clique size*/
    long unsigned int _decimation_move_var; /*number of variable that must be deciamted*/
    long unsigned int _updated_N; /*values of N updated with decimation moves*/

    struct _pred_sort_greater{/*predicator for Vertex sort*/
        
        bool operator()(Vertex *  &x,Vertex *  &y){
            
            return real(x->_message_i) > real(y->_message_i);
            
        }
        
    };
    
    struct _pred_sort_smaller_vertex{/*predicator for Vertex sort*/
        
        bool operator()(Vertex *  &x,Vertex *  &y){
            
            return x->_vertex < y->_vertex;
            
        }
        
    };
    
    template<class T>
    struct _complex_greater_pred {/*predicator for the complex sort*/
        
        bool operator()(const complex<T>  &x,const complex<T>  &y){
            
            return real(x) > real(y);
            
        }
        
    };

  
    /*private members*/
    long double _mean(vector<long double> &V){
        long double sum=0.;
        for (long unsigned int i=0; i<static_cast<long unsigned int>(V.size()); ++i) {
            sum+=static_cast<long double>(V[i]);
        }
        sum=sum/static_cast<long double>(V.size());
        return sum;
    }
    
    long double _std(vector<long double> &V, const long double &mean){
        long double sum=0.;
        for (long unsigned int i=0; i<static_cast<long unsigned int>(V.size()); ++i) {
            sum+=static_cast<long double>((V[i]-mean)*(V[i]-mean));
        }
        sum=sum/static_cast<long double>(V.size()-1);
        return sqrt(sum);
    }
    
    
    
    /******************************************/
    /******************************************/
    /******************************************/
    /*********** BP equations START ***********/
    /******************************************/
    /******************************************/
    /******************************************/
    
    /*START*/
    
    /*-3*/
    /*functions for computing part of BP equations*/
    
    inline long double _diff_logs(long double & _G_l_to_i, long unsigned int &i, long unsigned int &j){
        long double y=0;
        y=static_cast<long double>(log(1.+((1.+_cast_bool_to_int(i,j))*exp(_G_l_to_i)*_inv_sqrtN))-log(1.+(exp(_G_l_to_i)*_inv_sqrtN)));
        return y;
    }
    inline long double _diff_logs(long double &_G_l_to_i, unsigned int &i, unsigned int &j){
        long double y=0;
        y=static_cast<long double>(log(1.+((1.+_cast_bool_to_int(i,j))*exp(_G_l_to_i)*_inv_sqrtN))-log(1.+(exp(_G_l_to_i)*_inv_sqrtN)));
        return y;
    }
    inline long double _diff_logs(long double &_G_l_to_i, int &i, int &j){
        long double y=0;
        y=static_cast<long double>(log(1.+((1.+_cast_bool_to_int(i,j))*exp(_G_l_to_i)*_inv_sqrtN))-log(1.+(exp(_G_l_to_i)*_inv_sqrtN)));
        return y;
    }
    
    /*-2*/
    /*updating messages i*/
    
    void _sum_logs( vector<Vertex *> &V, long unsigned int &i){
        long double sum=0;
        long double _y=0;
        for (long unsigned int _l=0; _l<_N; ++_l) {
            if(_l!=i && V[_l]->_I_am_fixed){
                _y=real(V[_l]->_message_i_to_j[i]);
                sum+=_diff_logs(_y, _l, i);
            }
        }
        sum+=static_cast<long double>(log(_kappa));
        V[i]->_message_i=complex<long double>(real(V[i]->_message_i),sum);
    }
    
    /*-1*/
    
    /*updating messages*/
    bool _update_messages(vector<Vertex *> &V){
        long unsigned int counter_conv=0;
        long double _y=0;
        long double tempre=0, tempim=0;
        long double var=0.;
        bool flag=true;
        for (long unsigned int i=0; i<_N; ++i) {
            if(i!=V[i]->_vertex)exit(-1);
            if(V[i]->_I_am_fixed){
            _sum_logs(V, i);
            tempre=imag(V[i]->_message_i);
            tempim=real(V[i]->_message_i);
            V[i]->_message_i=complex<long double>(tempre,tempim);
                if((abs(tempre-tempim)<eps) && flag){
                    flag=true;
                    ++counter_conv;
                    
                }else{
                    flag=false;
                }

            }
        }

        /*updating messages i->j start*/
        for (long unsigned int i=0; i<_N; ++i) {
            if(V[i]->_I_am_fixed){
            for (long unsigned int j=0; j<_N; ++j) {
                 if(j!=i && V[j]->_I_am_fixed){
                    _y=real(V[j]->_message_i_to_j[i]);
                     var=real(V[i]->_message_i)-_diff_logs(_y, j, i);
                    V[i]->_message_i_to_j[j]=complex<long double>(real(V[i]->_message_i_to_j[j]),var);
                 }
            }
            }
        }
        for (long unsigned int i=0; i<_N; ++i) {
            if(V[i]->_I_am_fixed){
            for (long unsigned int j=0; j<_N; ++j) {
                if(V[j]->_I_am_fixed){
                tempre=imag(V[i]->_message_i_to_j[j]);
                tempim=real(V[i]->_message_i_to_j[j]);
                V[i]->_message_i_to_j[j]=complex<long double>(tempre,tempim);
                }
              
            }
        }
        }
        cout<<"counter_conv: "<<counter_conv<<" "<<flag<<endl;
        /*updating messages i->j end*/
        if(flag){
            return true;
        }else{
            return false;
        }
    }
    
    /*END*/
    
    /******************************************/
    /******************************************/
    /******************************************/
    /*********** BP equations END *************/
    /******************************************/
    /******************************************/
    /******************************************/
    
    void _update_Graph(){ /*Init values*/
        _frac_clique_found.reserve(NGraphs);
        MyRand.print_seed();
        order_Graph=_N;
        size_Graph=_M;
        size_HC=_HC;
        _inv_sqrtN=1./sqrt(_N);
        _counter=0;
        _ngbs.resize(order_Graph);
        _adjmatrix.resize(_N);
        _updated_N=_N;
        for (long unsigned int i=0; i<_N; ++i) {
            _adjmatrix[i].resize(_N);
        }
     //   _kappa=static_cast<long double>(1./sqrt(e));
        _kappa=static_cast<long double>((double)_HC/(double)sqrt(_N));
    }

    void _clean_ngbs(vector<Vertex *> &V){ /*clean vector ngbs*/
        _decimation_move_var=static_cast<long unsigned int>(prob_dec*((double)_N));
//       _kappa=static_cast<long double>(1./(double)sqrt(e));
         _kappa=static_cast<long double>((double)_HC/(double)sqrt(_N));
        _updated_N=_N;
        for (long unsigned int i=0; i<_N; ++i) {
            V[i]->_degree_vertex_i=0;
            V[i]->_message_i=complex<long double>(0.,0.);
            V[i]->_I_am_fixed=true;
            V[i]->_vertex=i;
            _ngbs[i].clear();
            _adjmatrix[i][i]=false;
        }
    }
};



int main(int argc, const char * argv[]) {
    vector<Vertex> V1;
    vector<Vertex *> V;
    /*Initialization vectors and  variables*/
    Graph Mygraph(argc,argv);
    V1.resize(Mygraph.order_Graph);
    V.resize(Mygraph.order_Graph);
    for (long unsigned int i=0; i<Mygraph.order_Graph; ++i) {
        V1[i]._vertex=i;
        V1[i]._message_i_to_j.resize(Mygraph.order_Graph);
        V[i]=&V1[i];
    }
    
    /*build and find the hidden clique*/
    clock_t begin, end;
    double elapsed_secsm=0., elapsed_secsm1=0.;
    
    for (int counterG=0; counterG<NGraphs; ++counterG) {
        Mygraph.build_graph(V); /*build a graph*/
        begin = clock();
        Mygraph._iterative_procedure(V); /*find the HC into a graph using AMP*/
        end = clock();
        elapsed_secsm+= double(end - begin) / CLOCKS_PER_SEC;
        begin = clock();
        Mygraph.SM_1_and_stop(V); /*find the HC into a graph  using SM1 early stopping*/
        end = clock();
        elapsed_secsm1+= double(end - begin) / CLOCKS_PER_SEC;
    }
    
    ostringstream n;
    ostringstream h;
    ostringstream hc;
    n<<Mygraph.order_Graph;
    h<<frac;
    hc<<Mygraph.give_me_HC();
    string s=directory;
    s+="mean";
    s+=n.str();
    s+="-f";
    s+=h.str();
    s+=".txt";
    ofstream outfile(s.c_str(), ios_base::app);
    outfile<<Mygraph;/*frac of sites of HC in the first HC positions*/
    s.clear();
    s=directory;
    s+="timing_AMP_VS_SM1.txt";
    ofstream outfiletime(s.c_str(), ios_base::app);
    outfiletime<<Mygraph.order_Graph<<" "<<Mygraph.size_HC<<" "<<elapsed_secsm<<" "<<elapsed_secsm1<<endl;
    s.clear();
    s=directory;
    s+="Readmetime.txt";
    ofstream readme(s.c_str());
    readme<<"N  HC  timeBP  timeSM1"<<endl;
    s.clear();
    s=directory;
    s+="Which_Graph";
    s+=n.str();
    s+="-hc";
    s+=hc.str();
    s+=".txt";
    ofstream out_which(s.c_str(), ios_base::app);
    for (unsigned long i=0; i<Mygraph._size_clique_hc_sites.size(); ++i) {
        out_which<<Mygraph.order_Graph<<" "<<Mygraph.give_me_HC()<<" "<<i<<" "<<real(Mygraph._size_clique_hc_sites[i])<<" "<<imag(Mygraph._size_clique_hc_sites[i])<<" "<<real(Mygraph._iteration_cleaning[i])<<" "<<imag(Mygraph._iteration_cleaning[i])<<" "<<real(Mygraph._size_clique_hc_sites_AMP[i])<<" "<<imag(Mygraph._size_clique_hc_sites_AMP[i])<<" "<<real(Mygraph._iteration_AMP[i])<<" "<<imag(Mygraph._iteration_AMP[i])<<endl;
    }
    
    return 0;
}
