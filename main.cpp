#include <iostream>
#include <cmath>
#include "minpack.h"
#include "gnuplot_i.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>
#include <windows.h>
#endif


using namespace std;


struct node{
    double price;
    double demand;
    double supply;
    node(){}
};


ostream& operator << (ostream & stream , const node & nd);
istream& operator >> (istream & stream , node & nd);
string print_model(double * a);
void set_plot_ranges(Gnuplot & plot);
void make_and_show_plot(Gnuplot & plot , double * supply_params , double * demand_params ,
                        pair<double , double> & equilibrium , char* path_to_exe);
void wait_for_key();
string inverse_model(double * params , double grant);


vector<node> t(12);

// y = a + b / 2 ^ (x / c)
void supply(int n, double x[], double y[], int & iflag){
    y[0] = y[1] = y[2] = 0;
    for(auto & e : t){
        double f = (x[0] + x[1] * pow(M_E , x[2] * e.price) - e.supply);
        y[0] += f;
        y[1] += f * pow(M_E , x[2] * e.price);
        y[2] += f * x[1] * e.price * pow(M_E , x[2] * e.price);
    }
    y[0] *= 2, y[1] *= 2, y[2] *= 2;
}

// y = a + b * e ^ cx
void demand(int n, double x[], double y[], int & iflag){
    y[0] = y[1] = y[2] = 0.0;
    for(auto & e : t){
        double f = (x[0] + x[1] * pow(M_E , x[2] * e.price) - e.demand);
        y[0] += f;
        y[1] += f * pow(M_E , x[2] * e.price);
        y[2] += f * x[1] * e.price * pow(M_E , x[2] * e.price);
    }
    y[0] *= 2, y[1] *= 2, y[2] *= 2;
}


double* estimate(void func(int , double[] , double[], int &), int n, vector<double> init_x){
    int cnt = 150 * n * n; //must be not less then (n * (3 * n + 13)) / 2
    double eps = 1e-7;
    double *x = new double[n];
    double *y = new double[n];
    double *wa = new double[cnt];

    for(int i = 0 ; i < n ; ++i) x[i] = init_x[i];

    hybrd1(func , n , x , y , eps , wa , cnt);

    delete [] y;
    delete [] wa;

    return x;
}

// y = a + b * e ^ cx
inline double supply_model(double * a , double x){
    return a[0] + a[1] * pow(M_E , a[2] * x);
}

// y = a + b * e ^ cx
inline double demand_model(double * a , double x){
    return a[0] + a[1] * pow(M_E , a[2] * x);
}


pair<double , double> equilibrium(double model1(double *, double),
                                  double model2(double *, double), double * a1 , double * a2){
    double lx = -1000 , rx = 1000;
    while(abs(rx - lx) >  1e-10){
        double midx = 0.5 * (lx + rx);
        double y1 = model1(a1 , lx);
        double y2 = model2(a2 , lx);
        double midy1 = model1(a1 , midx);
        double midy2 = model2(a2 , midx);

        if((y1 - y2 - 1e-10 > 0 && midy1 - midy2 - 1e-10 > 0) || (y1 - y2 + 1e-10 < 0 && midy1 - midy2 + 1e-10 < 0)){
            lx = midx;
        }else {
            rx = midx;
        }
    }
    return {lx , model1(a1 , lx)};
}

inline double avg_demand(){
    double sum = 0;
    for(auto & e : t) sum += e.demand;
    return sum / t.size();
}

inline double avg_supply(){
    double sum = 0;
    for(auto & e : t) sum += e.supply;
    return sum / t.size();
}

inline double avg_price(){
    double sum = 0;
    for(auto & e : t) sum += e.price;
    return sum / t.size();
}

double supply_segment_arc_elasticity(){
    double delta_q = t.back().supply - t.front().supply;
    double delta_p = t.back().price - t[0].price;
    return (delta_q / delta_p) / (avg_supply() / avg_price());
}

double demand_segment_arc_elasticity(){
    double delta_q = t.back().demand - t.front().demand;
    double delta_p = t.back().price - t.front().price;
    return (delta_q / delta_p) / (avg_demand() / avg_price());
}

string is_elastic(double k){
    assert(abs(k) != numeric_limits<double>::infinity());
    return abs(k) < 1 ? "inelastic" : "elastic";
}

double point_elasticity(double model(double *, double) , double * params , double price , double len = 0.001){
    double delta_p = len;
    double delta_q = model(params , price + len) - model(params , price);
    return (delta_q * price) / (delta_p * model(params , price));
}


const char* input_file = "D:\\Clion_projects\\MEEP_1\\input.txt";
string path_to_temp_files;

int main(int argc , char** argv) {
    assert(argc > 0);
    cout.setf(ios::fixed);
    cout.precision(5);
    freopen(input_file , "r" , stdin);

    for(auto & e : t) cin >> e;


/*
 * PLOT
 * path to 'gnuplot.exe' must be set at line 594 of gnupllot_i.hpp to Gnuplot::m_sGNUPlotPath
 */
    Gnuplot plot("lines");
    plot.set_title("MEEP Lab 1");


/*
 * ESTIMATING SUPPLY MODEL
 */
    double* sx = estimate(supply , 3 , {240 , -250 , -0.02});
    cout << endl << "Supply Model:" << endl;
    print_model(sx);


/*
 * ESTIMATING DEMAND MODEL
 */
    double* dx = estimate(demand, 3 , {10 , 200 , -0.2});
    cout << endl << "Demand model:" << endl;
    print_model(dx);


/*
 * FINDING AND EQUILIBRIUM POINT
 */
    cout << endl << "Equilibrium point" << endl;
    auto eq = equilibrium(supply_model, demand_model, sx, dx);
    cout << eq.first << " " << eq.second << endl << endl;


/*
 * CHECKING STABILITY OF AN EQUILIBRIUM
 */
    double sel = point_elasticity(supply_model , sx , eq.first);
    double dem = point_elasticity(demand_model , dx , eq.first);
    cout << "Supply elasticity at the equilibrium: " << sel << endl;
    cout << "Demand elasticity at the equilibrium: " << dem << endl;
    cout << (abs(dem) > abs(sel) ? "Stable" : "Unstable") << endl;


/*
 * ARC ELASTICITY
 */
    cout << endl << "Arc elasticity" << endl;
    double sup_el = supply_segment_arc_elasticity();
    double dem_el = demand_segment_arc_elasticity();
    cout << "Supply: " << sup_el << " - " << is_elastic(sup_el) << endl;
    cout << "Demand: " << dem_el << " - " << is_elastic(dem_el) << endl;


    make_and_show_plot(plot , sx , dx , eq , argv[0]);

    delete [] sx;
    delete [] dx;


    plot.remove_tmpfiles();
//    remove((path_to_temp_files + "\\supply.txt").c_str());
//    remove((path_to_temp_files + "\\demand.txt").c_str());


    return 0;
}


inline ostream& operator << (ostream & stream , const node & nd){
    stream << "Price " << nd.price << " , demand " << nd.demand
           << " , supply " << nd.supply;
}

inline istream& operator >> (istream & stream , node & nd){
    stream >> nd.price >> nd.demand >> nd.supply;
}

inline string print_model(double * a){
    cout << a[0] << (a[1] >= 0 ? "+" : "") << a[1] << "*e^(" << a[2] << "*x)";
    return to_string(a[0]) + (a[1] >= 0 ? "+" : "") + to_string(a[1]) + "*exp(1)**(" + to_string(a[2]) + "*x)";
}


//for an inverse model x = ln((y - a)/b) / c
void set_plot_ranges(Gnuplot & plot){
    double min_x = DBL_MAX , max_x = -DBL_MAX;
    double min_y = DBL_MAX , max_y = -DBL_MAX;
    for(auto & e : t){
        min_x = min(min_x , e.price);
        max_x = max(max_x , e.price);
        min_y = min(min_y , min(e.supply , e.demand));
        max_y = max(max_y , max(e.supply , e.demand));
    }
    plot.set_xrange(min_y - 1 , max_y + 1);
    plot.set_yrange(min_x - 1 , max_x + 1);
}


void make_and_show_plot(Gnuplot & plot , double * supply_params , double * demand_params ,
                        pair<double , double> & equilibrium , char* path_to_exe){
    set_plot_ranges(plot);
    plot.set_xlabel("Q");
    plot.set_ylabel("P");
    plot.set_grid();

    //supply and demand equations
    plot.plot_equation(inverse_model(supply_params , 0) , "Supply");
    plot.plot_equation(inverse_model(demand_params , 0) , "Demand");

    //print points to temp file
    string path_to_file(path_to_exe);
    path_to_file = path_to_file.substr(0 , path_to_file.find_last_of("/\\"));
    ofstream ofs(path_to_file + "\\supply.txt");
    for(auto & e : t){
        ofs << e.supply << "   " << e.price << endl;
    }
    ofs.close();

    ofs.open(path_to_file + "\\demand.txt");
    for(auto & e : t){
        ofs << e.demand << "   " << e.price << endl;
    }
    ofs.close();
    path_to_temp_files = path_to_file;

    //styles
    plot.cmd("set style line 1 lc rgb '#ac38dd' pt 7 ps 1.5 lt 1 lw 2\n");
    plot.cmd("set style line 3 lc rgb '#00ff90' pt 7 ps 1.5 lt 1 lw 2\n");
    plot.cmd("set style line 2 lc rgb '#ffc43a' pt 7 ps 1.5 lt 1 lw 2\n");

    //supply and demand points from table
    plot.cmd("replot '" + path_to_file + "\\supply.txt" + "' w p ls 1 title 'Supply points'\n");
    plot.cmd("replot '" + path_to_file + "\\demand.txt" + "' w p ls 3 title 'Demand points'\n");

    //equilibrium point
    plot.cmd("eqx=" + to_string(equilibrium.second) + "\n");
    plot.cmd("eqy=" + to_string(equilibrium.first) + "\n");
    plot.cmd("replot '+' using (eqx):(eqy) title \"Equilibrium\"  w p ls 2\n");

    //grant
    for(int log2g = 0 ; log2g < 4 ; ++log2g){
        plot.plot_equation(inverse_model(supply_params , 1 << log2g)
                , "Supply with grant of size " + to_string(1 << log2g));
    }

    plot.showonscreen();
    wait_for_key();
}

void wait_for_key(){
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    cout << endl << "Press ENTER to continue..." << endl;
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    cin.clear();
    cin.ignore(cin.rdbuf()->in_avail());
    cin.get();
#endif
    return;
}

string inverse_model(double * params , double grant){
    string s = "((log((x-" + to_string(params[0]) + ")/" + to_string(params[1]) +
            "))/" + to_string(params[2]) + ") - " + to_string(grant);
    return s;
}
