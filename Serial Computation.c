# include <stdio.h>
# include <time.h> 
# include <math.h> 
# include <stdlib.h>

int main(int argc, char * argv[]);
void compute(int np, int nd, double pos[], double vel[], double mass, double f[], double * pot, double * kin);
double dist(int nd, double r1[], double r2[], double dr[]);
void initialize(int np, int nd, double box[], int * seed, double pos[], double vel[], double acc[]);
double r8_uniform_01(int * seed);
void timestamp(void);
void update(int np, int nd, double pos[], double vel[], double f[], double acc[], double mass, double dt);

int main(int argc, char * argv[])
{
    clock_t start_time, end_time;
    double * acc;
    double * box;
    double dt = 0.0001;
    double e0;
    double * force;
    int i;
    int id;
    double kinetic;
    double mass = 1.0;
    int nd = 3;
    int np = 1000;
    double * pos;
    double potential;
    int proc_num;
    int seed = 123456789;
    int step;
    int step_num = 400;
    int step_print;
    int step_print_index;
    int step_print_num;
    double * vel;
    double wtime;
    double c=299792458;
    start_time = clock();
    timestamp();
    acc = (double *)malloc(nd * np * sizeof (double));
    box = (double *)malloc(nd * sizeof(double));
    force = (double *)malloc(nd * np * sizeof(double));
    pos = (double *)malloc(nd * np * sizeof(double));
    vel = (double *)malloc(nd * np * sizeof(double));
    printf("\n");
    
    printf("-------------------------------------Molecular dynamics program------------------------------------\n");
    printf("\n");
    printf("Number of particles in the simulation: %d\n", np);
    printf("\n");
    printf("Number of time steps: %d\n", step_num);
    printf("\n");
    printf("Size of each time step: %f\n", dt);
    printf("\n");
    for (i = 0; i < nd; i++) {
        box[i] = 10.0;
    }
    printf("\n");
    //printf("Initializing positions, velocities, and accelerations.\n");
    initialize(np, nd, box, & seed, pos, vel, acc);
    printf("\n");
    //printf("Computing initial forces and energies.\n");
    compute(np, nd, pos, vel, mass, force, & potential, & kinetic);
    printf("POS:%f \nVelocity:%f \nMass:%f \nForce:%f \n", *pos, *vel, mass, *force);
    e0 = potential + kinetic;
    printf("\n");
   
    printf("\n");
    printf("Step\tPotential Energy\tKinetic Energy\tRelative Energy Error\tKinetic Energy Velocity\n");
    printf("----\t----------------\t--------------\t----------------------\t---------------------\n");
    printf("\n");
    step_print = 0;
    step_print_index = 0;
    step_print_num = 10;
    step = 0;
    printf("%4d\t%16f\t%13f\t%21e\t%22e\n", step, potential, kinetic, (potential + kinetic - e0) / e0, (sqrt(2*kinetic)/mass));
    step_print_index = step_print_index + 1;
    step_print = (step_print_index * step_num) / step_print_num;
    for (step = 1; step <= step_num; step++) {
        compute(np, nd, pos, vel, mass, force, & potential, & kinetic);
        if (step == step_print) {
           
            printf("%4d\t%16f\t%13f\t%21e\t%22e\n", step, potential, kinetic, (potential + kinetic - e0) / e0, (sqrt(2*kinetic)/mass));
            step_print_index = step_print_index + 1;
            step_print = (step_print_index * step_num) / step_print_num;
        }
        update(np, nd, pos, vel, force, acc, mass, dt);
    }
    printf("\n\n");
    printf("Step\tThermal Energy\t\tMomentum\tElectromagnetic Energy\tRelativistic Energy\n");
    printf("----\t----------------\t--------------\t----------------------\t---------------------\n");
    printf("\n");
    step_print = 0;
    step_print_index = 0;
    step_print_num = 10;
    step = 0;
    double velocity = *vel;
    printf("%4d\t%16f\t%13f\t%21e\t%22e\n", step, kinetic, potential, (sqrt(2*kinetic*mass)),(mass * c * c / sqrt(1 - ((velocity) * (velocity) / (c * c)))));
    step_print_index = step_print_index + 1;
    step_print = (step_print_index * step_num) / step_print_num;
    for (step = 1; step <= step_num; step++) {
        compute(np, nd, pos, vel, mass, force, & potential, & kinetic);
        velocity = *vel;
        if (step == step_print) {
        printf("%4d\t%16f\t%13f\t%21e\t%22e\n", step,kinetic,potential, (sqrt(2*kinetic*mass)),(mass * c * c / sqrt(1 - ((velocity) * (velocity) / (c * c)))));
            step_print_index = step_print_index + 1;
            step_print = (step_print_index * step_num) / step_print_num;
        }
        update(np, nd, pos, vel, force, acc, mass, dt);
    }
    printf("\n");
    end_time = clock();
    wtime = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Elapsed time for main computation:\n");
    printf("%f seconds.\n", wtime);
    free(acc);
    free(box);
    free(force);
    free(pos);
    free(vel);
    printf("\n");
    //MD_OPENMP
    //Normal end of execution
    printf("\n");
    timestamp();
    return 0;
}
/***/
void compute(int np, int nd, double pos[], double vel[], double mass, double f[], double * pot, double * kin)
{
    double d;
    double d2;
    int i;
    int j;
    int k;
    double ke;
    double pe;
    double PI2 = 3.141592653589793 / 2.0;
    double rij[3];
    pe = 0.0;
    ke = 0.0;
    for (k = 0; k < np; k++) {
        for (i = 0; i < nd; i++) {
            f[i + k * nd] = 0.0;
        }
        for (j = 0; j < np; j++) {
            if (k != j) {
                d = dist(nd, pos + k * nd, pos + j * nd, rij);
                if (d < PI2) {
                    d2 = d;
                } else {
                    d2 = PI2;
                }
                pe = pe + 0.5 * pow(sin(d2), 2);
                for (i = 0; i < nd; i++) {
                    f[i + k * nd] = f[i + k * nd] - rij[i] * sin(2.0 * d2) / d;
                }
            }
        }
        for (i = 0; i < nd; i++) {
            ke = ke + vel[i + k * nd] * vel[i + k * nd];
        }
    }
    ke = ke * 0.5 * mass;
    * pot = pe;
    * kin = ke;
    return;
}
/****/

double dist(int nd, double r1[], double r2[], double dr[])
{
    double d;
    int i;
    d = 0.0;
    for (i = 0; i < nd; i++) {
        dr[i] = r1[i] - r2[i];
        d = d + dr[i] * dr[i];
    }
    d = sqrt(d);
    return d;
}
/****/

void initialize(int np, int nd, double box[], int * seed, double pos[], double vel[], double acc[])
{
    int i;
    int j;
    for (i = 0; i < nd; i++) {
        for (j = 0; j < np; j++) {
            pos[i + j * nd] = box[i] * r8_uniform_01(seed);
        }
    }
    for (j = 0; j < np; j++) {
        for (i = 0; i < nd; i++) {
            vel[i + j * nd] = 0.0;
        }
    }
    for (j = 0; j < np; j++) {
        for (i = 0; i < nd; i++) {
            acc[i + j * nd] = 0.0;
        }
    }
    return;
}
/****/

double r8_uniform_01(int * seed)
{
    int k;
    double r;
    k = * seed / 127773;
    * seed = 16807 * ( * seed - k * 127773) - k * 2836;
    if ( * seed < 0) {
        * seed = * seed + 2147483647;
    }
    r = (double)( * seed) * 4.656612875E-10;
    return r;
}
/****/

void timestamp(void)
{
    # define TIME_SIZE 40
    static char time_buffer[TIME_SIZE];
    const struct tm * tm;
    size_t len;
    time_t now;
    now = time(NULL);
    tm = localtime( & now);
    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
    printf("%s\n", time_buffer);
    return;
    # undef TIME_SIZE
}
/***/

void update(int np, int nd, double pos[], double vel[], double f[], double acc[], double mass, double dt)
{
    int i;
    int j;
    double rmass;
    rmass = 1.0 / mass;
    for (j = 0; j < np; j++) {
        for (i = 0; i < nd; i++) {
            pos[i + j * nd] = pos[i + j * nd] + vel[i + j * nd] * dt + 0.5 * acc[i + j * nd] * dt * dt;
            vel[i + j * nd] = vel[i + j * nd] + 0.5 * dt * (f[i + j * nd] * rmass + acc[i + j * nd]);
            acc[i + j * nd] = f[i + j * nd] * rmass;
        }
    }
    return;
}