#include <climits>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
#include <queue>
#include <climits>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/basic_graph_alg.h>
#include <LEDA/graph/min_cut.h>
#include <LEDA/graph/max_flow.h>
#include <set>
using namespace std;
using namespace leda;

struct Gate{
    std::string name;
    leda::string lName;
    bool isPI;
    bool isPO;
    std::vector<std::string> table;
    std::vector<std::string> inputs;
    node gNode;
    int label;
    
    Gate(){};
    Gate(std::string name, leda::string lName, bool isPI, bool isPO,
            std::vector<std::string> tb, std::vector<std::string> inps, node gNode, int label) : 
        name(name), lName(lName), isPI(isPI), isPO(isPO), table(tb), inputs(inps), gNode(gNode), label(label){}
};

using GateMap = std::map<std::string, Gate>;
using NameMap = std::map<node, std::string>;
using EdgeWeight = std::map<edge, int>;
using Cluster = std::map<std::string, std::vector<std::string>>;

std::string graphName;
std::vector<Gate> primeInput;
std::vector<Gate> primeOutput;
graph bigGraph;
graph outputGraph;
GateMap outGM;
NameMap outNM;
std::vector<edge> edges;
GateMap gateMap;
NameMap nameMap;
Cluster cluster;

int K_CUT;


void readBLIF(char * fileName);
void decomposeMultiInputGate();
int updateLabel(GateMap & gateMap, NameMap & nameMap, graph & gra);
void printGate(GateMap & gateMap);
void printGate(GateMap & gateMap, NameMap & nameMap, graph & g);
void printEdge();
void printEdge(EdgeWeight & ew, GateMap & gateMap);
void printPrime(std::vector<Gate> & gl, std::string p);
void writeBLIF(char * fname);
void makeSubGraph(graph & sg, node root, GateMap & subGateMap, NameMap & subNameMap);
void flowMapPhase1();
int mergeMaxLabelNode(graph & sg, node root, GateMap & subGateMap, NameMap & subNameMap);
node addSource(graph & sg, GateMap * subGateMap, NameMap & subNameMap);
void phase1TransferGraph(graph & sg, node source, node root, GateMap & subGateMap, NameMap & subNameMap);
bool findMinCut(graph & sg, node source, node root, GateMap & subGateMap, NameMap & subNameMap);
void flowMapPhase2();

int main(int argc, char ** argv){
    
    if(argc != 3){
        return -1;
    }
    //cout << gateMap.size() << endl;
    readBLIF(argv[1]);
    decomposeMultiInputGate();
    flowMapPhase1();
    flowMapPhase2();
    printGate(outGM);
    outputGraph.print();
    //printGate(gateMap);
    int label = updateLabel(outGM, outNM, outputGraph);
    //printGate(outGM);
    cout << "after mapping label: " << label << endl;

    //for(auto & iter : cluster){
        //cout << "name : " << iter.first << endl;
        //cout << "cluster: ";
        //for(auto & n : iter.second){
            //cout << n << " ";
        //}
        //cout << "\n========================\n";
    //}
    //printGate(gateMap);
    //bigGraph.print();
    //printPrime(primeInput, "input");
    //writeBLIF(argv[2]);
    //printPrime(primeOutput);

    //node_array<int> ord(bigGraph);
    //TOPSORT(bigGraph, ord);
    //node v;
    //forall_nodes(v, bigGraph){
        //cout << gateMap[nameMap[v]].name << " ";
        //cout << ord[v] << endl;
    //}
    //for(auto & iter : gateMap){
        //cout << "name: " << iter.second.name << endl;
        //cout << "label: " << iter.second.label << endl;
    //}

    //bigGraph.print();
    //for(int i = 0;i < primeInput.size(); ++i){
        //cout << nameMap[primeInput[i].gNode] << endl;
    //}
    

    return 0;
}

void flowMapPhase2(){
    std::vector<Gate> mappingList = primeOutput;
    while(!mappingList.empty()){
        int length = mappingList.size();
        std::set<std::string> inputList;
        Gate & g = gateMap[mappingList[0].name];
        auto & c = cluster[g.name];
        for(auto & cName : c){
            //cout << cName << endl;
            Gate & gg = gateMap[cName];
            leda::list<edge> elist = bigGraph.in_edges(gg.gNode);
            while(!elist.empty()){
                edge e = elist.Pop();
                node n = e->terminal(0);
                if(gateMap[nameMap[n]].label < g.label){
                    if(!gateMap[nameMap[n]].isPI){
                        bool add = true;
                        for(int i = 0; i < length; ++i){
                            Gate & og = mappingList[i];
                            if(og.name == nameMap[n]){
                                add = false;
                                break;
                            }
                        }
                        if(add){
                            mappingList.push_back(gateMap[nameMap[n]]);
                        }
                        inputList.insert(nameMap[n] + "LUT");
                    }
                    else{
                        Gate pi = gateMap[nameMap[n]];
                        pi.name = pi.name + "LUT";
                        if(outGM.find(nameMap[n] + "LUT") == outGM.end()){
                            pi.gNode = outputGraph.new_node();
                            outGM[pi.name] = pi;
                            outNM[pi.gNode] = pi.name;
                        }
                        inputList.insert(pi.name);
                    }
                }
            }
        }
        Gate lut(g.name + "LUT", g.lName + "LUT", g.isPI, g.isPO, c, {}, outputGraph.new_node(), -1);
        std::copy(inputList.begin(), inputList.end(), std::back_inserter(lut.inputs));
        outGM[lut.name] = lut;
        outNM[lut.gNode] = lut.name;
        mappingList.erase(mappingList.begin());
    }
   
    node v;
    forall_nodes(v, outputGraph){
        Gate & g = outGM[outNM[v]];
        auto & ns = g.inputs;
        for(auto & n : ns){
            Gate & gg = outGM[n];
            outputGraph.new_edge(gg.gNode, g.gNode);
        }
    }
}

bool findMinCut(graph & sg, node source, node root, GateMap & subGateMap, NameMap & subNameMap){
    root = subGateMap[nameMap[root]].gNode;
    leda::edge_array<int> weight(sg);
    //add weight
    leda::list<edge> edgeList = sg.all_edges();
    while(!edgeList.empty()){
        //cout << "is in\n";
        edge e = edgeList.Pop();
        node from = e->terminal(0);
        node to = e->terminal(1);
        Gate & fg = subGateMap[subNameMap[from]];
        Gate & tg = subGateMap[subNameMap[to]];
        //cout << fg.name << " " << tg.name << endl;
        if(fg.name == tg.name + "p"){
            //cout << "clone gate: " << fg.name << endl;
            weight[e] = 1;
        }
        else{
            weight[e] = 1000;
        }
    }
    K_CUT = 3;
    leda::edge_array<int> flow;
    //edge ee;
    //forall_edges(ee, sg){
        //cout << weight[ee] << endl;
    //}
    int cut_value = MAX_FLOW(sg, source, root, weight, flow);
    edge e;
    //if(subNameMap[root] == "y1"){
        //cout << "y1 cut value: " << cut_value << endl;
    //}
    return cut_value <= K_CUT;
}

int mergeMaxLabelNode(graph & sg, node root, GateMap & subGateMap, NameMap & subNameMap){
    root = subGateMap[nameMap[root]].gNode;
    std::queue<node> nodeQueue;
    nodeQueue.push(root);
    int maxLabel = 0;
    leda::list<edge> el = sg.in_edges(root);
    GateMap tmpGateMap;
    while(!el.empty()){
        edge e = el.Pop();
        node n = e->terminal(0);
        maxLabel = subGateMap[subNameMap[n]].label > maxLabel ? subGateMap[subNameMap[n]].label : maxLabel;
    }
    while(!nodeQueue.empty()){
        node n = nodeQueue.front();
        nodeQueue.pop();
        leda::list<edge> iedgeList = sg.in_edges(n);
        while(!iedgeList.empty()){
            edge e = iedgeList.Pop();
            node nn = e->terminal(0);
            if(tmpGateMap.find(subNameMap[nn]) == tmpGateMap.end()){
                nodeQueue.push(nn);
                tmpGateMap[subNameMap[nn]] = subGateMap[subNameMap[nn]];
            }
        }
        leda::list<edge> edgeList = sg.out_edges(n);
        Gate & g = subGateMap[subNameMap[n]];
        while(!edgeList.empty()){
            edge e = edgeList.Pop();
            node nn = e->terminal(1);
            if(subGateMap[subNameMap[nn]].label == maxLabel){
                //if(subNameMap[root] == "y1"){
                    //cout << "label y1: " << subNameMap[nn] << endl;
                //}
                sg.new_edge(n, root);
                break;
            }
        }
    }

    //delete node
    nodeQueue.push(root);
    while(!nodeQueue.empty()){
        node n = nodeQueue.front();
        nodeQueue.pop();
        leda::list<edge> edgeList = sg.in_edges(n);
        while(!edgeList.empty()){
            edge e = edgeList.Pop();
            node nn = e->terminal(0);
            nodeQueue.push(nn);
        }
        Gate & g = subGateMap[subNameMap[n]];
        if(g.label == maxLabel){
            //if(subNameMap[root] == "y1"){
                //cout << "label y1: " << g.name << endl;
            //}
            cluster[subNameMap[root]].push_back(g.name);
            sg.del_node(n);
            subNameMap.erase(subNameMap.find(g.gNode));
            subGateMap.erase(subGateMap.find(g.name));
        }
    }
    return maxLabel;
}

node addSource(graph & sg, GateMap & subGateMap, NameMap & subNameMap){
    node source = sg.new_node();
    for(auto g : subGateMap){
        if(g.second.isPI){
            sg.new_edge(source, g.second.gNode);
        }
    }
    Gate sourceGate("source", "source", false, false, {}, {}, source, -1);
    subGateMap["source"] = sourceGate;
    subNameMap[source] = "source";
    return source;
}

void phase1TransferGraph(graph & sg, node source, node root, GateMap & subGateMap, NameMap & subNameMap){
    root = subGateMap[nameMap[root]].gNode;
    std::queue<node> nodeQueue;
    nodeQueue.push(root);
    GateMap tmpNodeList;
    while(!nodeQueue.empty()){
        node n = nodeQueue.front();
        nodeQueue.pop();
        leda::list<edge> el = sg.in_edges(n);
        while(!el.empty()){
            edge e = el.Pop();
            node nn = e->terminal(0);
            if(tmpNodeList.find(subNameMap[nn]) == tmpNodeList.end()){
                nodeQueue.push(nn);
                tmpNodeList[subNameMap[nn]] = subGateMap[subNameMap[nn]];
            }
        }
        if(n != source && n != root){
            //add new node to self and set edge weight;
            Gate & ng = subGateMap[subNameMap[n]];
            //cout << "name: " << ng.name << endl;;
            node np = sg.new_node();
            Gate ngp(ng.name+"p", ng.lName + "p", ng.isPI, ng.isPO, {}, {}, np, ng.label);
            subGateMap[ngp.name] = ngp;
            subNameMap[ngp.gNode] = ngp.name;
            leda::list<edge> inEdge = sg.in_edges(n);
            while(!inEdge.empty()){
                edge e = inEdge.Pop();
                node nn = e->terminal(0);
                sg.new_edge(nn, np);
            }
            inEdge = sg.in_edges(n);
            sg.del_edges(inEdge);
            sg.new_edge(np, n);
        }
    }
}

void flowMapPhase1(){
    node_array<int> ord(bigGraph);
    TOPSORT(bigGraph, ord);
    node v;
    std::vector<node> topoOrderNode(gateMap.size());
    forall_nodes(v, bigGraph){
        topoOrderNode[ord[v] - 1] = v;
    }
    for(auto & n : topoOrderNode){
        if(!gateMap[nameMap[n]].isPI){
            graph sg;
            GateMap subGM;
            NameMap subNM;
            makeSubGraph(sg, n, subGM, subNM);
            node source = addSource(sg, subGM, subNM);
            int p = mergeMaxLabelNode(sg, n, subGM, subNM);
            //if(nameMap[n] == "y1"){
                //printGate(subGM);
            //}
            phase1TransferGraph(sg, source, n, subGM, subNM);
            if(findMinCut(sg, source, n, subGM, subNM)){
                gateMap[nameMap[n]].label = p;
                cluster[nameMap[n]].push_back(nameMap[n]);
            }            
            else{
                gateMap[nameMap[n]].label = p + 1;
                cluster[nameMap[n]] = {nameMap[n]};
            }
            //if(nameMap[n] == "y1"){
                //printGate(subGM, subNM, sg);
                //sg.print();
            //}
        }
    }

}

void makeSubGraph(graph & sg, node root, GateMap & subGateMap, NameMap & subNameMap){
    std::queue<node> nodeQueue;
    nodeQueue.push(root); // from big graph
    while(!nodeQueue.empty()){
        node n = nodeQueue.front();
        nodeQueue.pop();
        Gate g;
        if(subGateMap.empty()){
            //handle root;
            g = gateMap[nameMap[n]];
            g.gNode = sg.new_node();
            g.isPO = true;
            subGateMap[g.name] = g;
            subNameMap[g.gNode] = g.name;
        }
        else{
            g = subGateMap[gateMap[nameMap[n]].name];
        }
        leda::list<edge> el = bigGraph.in_edges(n);
        while(!el.empty()){
            edge ee = el.Pop();
            node nn = ee->terminal(0);
            if(subGateMap.find(nameMap[nn]) != subGateMap.end()){
                sg.new_edge(subGateMap[nameMap[nn]].gNode, g.gNode);
            }
            else{
                nodeQueue.push(nn);
                Gate gg = gateMap[nameMap[nn]];
                gg.gNode = sg.new_node();
                subGateMap[gg.name] = gg;
                subNameMap[gg.gNode] = gg.name;
                sg.new_edge(gg.gNode, g.gNode);
            }
        }
    }
}


void writeBLIF(char * fname){
    node_array<int> ord(bigGraph);
    TOPSORT(bigGraph, ord);
    node v;
    std::vector<node> topoOrderNode(gateMap.size());
    forall_nodes(v, bigGraph){
        topoOrderNode[ord[v] - 1] = v;
    }
    ofstream outputfile(fname);
    outputfile << ".model " +  graphName << endl;
    outputfile << ".inputs "; 
    for(auto & g : primeInput){
        outputfile << g.name << " ";
    }
    outputfile << "\n.outputs ";
    for(auto & g : primeOutput){
        outputfile << g.name << " ";
    }
    outputfile << endl;
    for(auto & n : topoOrderNode){
        Gate & g = gateMap[nameMap[n]];
        if(!g.isPI){
            outputfile << ".names ";
            for(auto & n : g.inputs){
                outputfile << n << " ";
            }
            outputfile << g.name << endl;
            for(auto & t : g.table){
                outputfile << t << endl;
            }
        }
    }
    outputfile << ".end";
}

void printPrime(std::vector<Gate> & gl, std::string p){
    cout << "prime " + p + ":\n";
    for(auto & g : gl){
        cout << g.name << endl;
    }
    cout << "===============\n";
}

void printEdge(EdgeWeight & ew, GateMap & gateMap){
    for(auto & ed : ew){
        edge e = ed.first;
        cout << "from: " << gateMap[nameMap[(e->terminal(0))]].name << " to " << gateMap[nameMap[(e->terminal(1))]].name << endl;
    }
}

void printEdge(){
    for(int i = 0;i < edges.size(); ++i){
        edge e = edges[i];
        cout << "from: " << gateMap[nameMap[(e->terminal(0))]].name << " to " << gateMap[nameMap[(e->terminal(1))]].name << endl;
    }
}

void printGate(GateMap & gateMap, NameMap & nameMap, graph & g){
    node n;
    cout << "********************\n";
    forall_nodes(n, g){
        Gate & gg = gateMap[nameMap[n]];
        cout << "Gate info:\n";
        cout << "name: " << gg.name << endl;
        g.print_node(n);
        cout << endl;
        cout << "==============\n";
    }
    cout << "********************\n";
}

void printGate(GateMap & gateMap){
    
    cout << "********************\n";
    for(auto & g : gateMap){
        Gate & gg = g.second;
        cout << "Gate info:\n";
        cout << "name: " << gg.name << endl;;
        cout << "inputs: ";
        for_each(gg.inputs.begin(), gg.inputs.end(), [](std::string & n){
                    cout << n << " ";
                });
        cout << endl;
        cout << "table:\n";
        for_each(gg.table.begin(), gg.table.end(), [](std::string & n){
                    cout << n << "\n";
                });
        cout << "label: " << gg.label << endl;
        cout << "==============\n";
    }
    cout << "********************\n";
}

void readBLIF(char * fileName){
    ifstream inputBLIF(fileName);
    std::string line;
    Gate * curGate = nullptr;
    while(getline(inputBLIF, line)){
        std::stringstream ss(line);
        std::string mode;
        ss >> mode;
        if(mode == ".model"){
            ss >> graphName;
        }
        else if(mode == ".end"){
            break;
        }
        else if(mode == ".inputs"){
            std::string name;
            while(ss >> name){
                if(name == "\\"){
                    getline(inputBLIF, line);
                    ss.str(line);
                    ss.clear();
                }
                else{
                    Gate gate;
                    gate.name = name;
                    gate.lName = name.data();
                    gate.isPI = true;
                    gate.isPO = false;
                    gate.gNode = bigGraph.new_node();
                    gate.label = 0;
                    gateMap[name] = gate;
                    nameMap[gate.gNode] = name;
                    //cout << bigGraph.get_node_entry_string(gate.gNode).c_str() << endl;
                    //cout << name << endl;
                    primeInput.push_back(gate);
                }
            }
        }
        else if(mode == ".outputs"){
            std::string name;
            while(ss >> name){
                if(name == "\\"){
                    getline(inputBLIF, line);
                    ss.str(line);
                    ss.clear();
                }
                else{
                    Gate gate;
                    gate.name = name;
                    gate.lName = name.data();
                    gate.isPI = false;
                    gate.isPO = true;
                    gate.gNode = bigGraph.new_node();
                    gateMap[name] = gate;
                    nameMap[gate.gNode] = name;
                    primeOutput.push_back(gate);
                }
            }
        }
        else if(mode == ".names"){
            std::string input;
            std::vector<std::string> inputs;
            while(ss >> input){
                if(input == "\\"){
                    cout << "fuck damn\n";
                }
                inputs.push_back(input);
            }
            std::string name;
            name = inputs.back();
            inputs.pop_back();
            
            if(gateMap.find(name) != gateMap.end()){
                //primeout
                curGate = &gateMap[name];
                curGate->inputs = inputs;
            }
            else{
                Gate g;
                g.name = name;
                g.lName = name.c_str();
                g.isPI = false;
                g.isPO = false;
                g.gNode = bigGraph.new_node();
                g.inputs = inputs;
                gateMap[name] = g;
                nameMap[g.gNode] = name;
                curGate = &gateMap[name];
            }
            int labelMax = 0;
            for(auto & n : inputs){
                labelMax = gateMap[n].label > labelMax ? gateMap[n].label : labelMax;
                edge e = bigGraph.new_edge(gateMap[n].gNode, curGate->gNode);
                edges.push_back(e);
            }
            curGate->label = labelMax + 1;
        }
        else{
            if(curGate != nullptr){
                curGate->table.push_back(line);
            }
            else{
                cout << "need gate to have truth table\n";
                return;
            }
        }
        //remain edge weight depth
    }
    curGate = nullptr;
    inputBLIF.close();
}

void decomposeMultiInputGate(){
    
    node_array<int> ord(bigGraph);
    TOPSORT(bigGraph, ord);
    node v;
    std::vector<node> topoOrderNode(gateMap.size());
    forall_nodes(v, bigGraph){
        topoOrderNode[ord[v] - 1] = v;
    }

    auto f = [&](node & a, node & b){
        return gateMap[nameMap[a]].label < gateMap[nameMap[b]].label;
    };

    for(auto & i : topoOrderNode){
        auto es = bigGraph.in_edges(i);
        std::list<node> sortNodes;
        if(es.size() > 2){
            //cout << "muti input node: " << gateMap[nameMap[i]].name << endl;
            //cout << "node: " << endl;
            Gate & og = gateMap[nameMap[i]];
            bigGraph.del_edges(es);
            while(es.size() != 0){
                edge e = es.Pop();
                auto iter = std::find(edges.begin(), edges.end(), e);
                edges.erase(iter);
                sortNodes.push_back(e->terminal(0));
            }
            sortNodes.sort(f);
            int count = 1;
            while(!sortNodes.empty()){
                std::string newName = og.name + "sub" + std::to_string(count);
                node n1 = sortNodes.front();
                sortNodes.pop_front();
                node n2 = sortNodes.front();
                sortNodes.pop_front();
                std::string name1 = gateMap[nameMap[n1]].name;
                std::string name2 = gateMap[nameMap[n2]].name;
                if(og.table.size() > 1){
                    //or gate
                    std::string in1;
                    std::string in2;
                    if(sortNodes.empty()){
                        og.inputs = {name1, name2};
                        og.table = {"1- 1", "-1 1"};
                        og.label = std::max(gateMap[nameMap[n1]].label, gateMap[nameMap[n2]].label) + 1;
                        edge e1 = bigGraph.new_edge(n1, og.gNode);
                        edge e2 = bigGraph.new_edge(n2, og.gNode);
                        edges.push_back(e1);
                        edges.push_back(e2);
                    }
                    else{
                        int i;
                        for(i = 0; i < og.inputs.size();++i){
                            if(name1 == og.inputs[i]){
                                in1 = in1 + og.table[i][i] + "- 1";
                                break;
                            }
                        }
                        if(i >= og.inputs.size()){
                            in1 += "1- 1";
                        }
                        for(i = 0; i < og.inputs.size();++i){
                            if(name2 == og.inputs[i]){
                                in2 = in2 + "-" + og.table[i][i] + " 1";
                                break;
                            }
                        }
                        if(i >= og.inputs.size()){
                            in2 += "-1 1";
                        }
                        Gate newg(newName, newName.c_str(), false, false, {in1, in2}, {name1, name2}, bigGraph.new_node(), -1);
                        newg.label = std::max(gateMap[nameMap[n1]].label, gateMap[nameMap[n2]].label) + 1;
                        edge e1 = bigGraph.new_edge(n1, newg.gNode);
                        edge e2 = bigGraph.new_edge(n2, newg.gNode);
                        edges.push_back(e1);
                        edges.push_back(e2);
                        nameMap[newg.gNode] = newName;
                        gateMap[newName] = newg;
                        sortNodes.push_back(newg.gNode);
                    }
                }
                else{
                    //and gate
                    std::string in;
                    if(sortNodes.empty()){
                        og.inputs = {name1, name2};
                        og.table = {"11 1"};
                        og.label = std::max(gateMap[nameMap[n1]].label, gateMap[nameMap[n2]].label) + 1;
                        edge e1 = bigGraph.new_edge(n1, og.gNode);
                        edge e2 = bigGraph.new_edge(n2, og.gNode);
                        edges.push_back(e1);
                        edges.push_back(e2);
                    }
                    else{
                        int i;
                        for(i = 0;i < og.inputs.size();++i){
                            if(name1 == og.inputs[i]){
                                in += og.table[0][i];
                                break;
                            }
                        }
                        if(i >= og.inputs.size()){
                            in += "1";
                        }
                        for(i = 0;i < og.inputs.size();++i){
                            if(name2 == og.inputs[i]){
                                in += og.table[0][i];
                                break;
                            }
                        }
                        if(i >= og.inputs.size()){
                            in += "1";
                        }
                        in += " 1";
                        Gate newg(newName, newName.c_str(), false, false, {in}, {name1, name2}, bigGraph.new_node(), -1);
                        newg.label = std::max(gateMap[nameMap[n1]].label, gateMap[nameMap[n2]].label) + 1;
                        edge e1 = bigGraph.new_edge(n1, newg.gNode);
                        edge e2 = bigGraph.new_edge(n2, newg.gNode);
                        edges.push_back(e1);
                        edges.push_back(e2);
                        nameMap[newg.gNode] = newName;
                        gateMap[newName] = newg;
                        sortNodes.push_back(newg.gNode);
                    }
                }
                count += 1;
                sortNodes.sort(f);
            }
            //for(auto & n : sortNodes){
                //cout << gateMap[nameMap[n]].name << " " << gateMap[nameMap[n]].label << endl;
            //}
            //cout << "===========\n";

            //cout << "mutiple input gate should be decompose: " << gateMap[nameMap[i]].name << endl;
        }
    }
    int label = updateLabel(gateMap, nameMap, bigGraph);
    cout << "before mapping label: " << label << endl;
}

int updateLabel(GateMap & gateMap, NameMap & nameMap, graph & gra){
    node_array<int> ord(gra);
    TOPSORT(gra, ord);
    node v;
    std::vector<node> topoOrderNode(gateMap.size());
    forall_nodes(v, gra){
        topoOrderNode[ord[v] - 1] = v;
    }
    int max = 0;
    for(auto & n : topoOrderNode){
        Gate & g = gateMap[nameMap[n]];
        if(!g.isPI){
            int maxl = 0;
            for(auto & name : g.inputs){
                if(gateMap[name].label > maxl){
                    maxl = gateMap[name].label;
                }
            }
            g.label = maxl + 1;
            max = g.label;
        }
    }
    return max;
}
