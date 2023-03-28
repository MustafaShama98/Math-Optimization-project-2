#include <stdio.h>
#include <math.h>
#define NMAX 100
#define M_PI 3.141592

typedef double (*FUN_PTR)(double[]);
double h0(double x[]);
double h00(double x[]);
double h1(double x[]);
typedef void (*GRAD_FUN_PTR)(double grad[], double x[]);
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);

FUN_PTR objective_function;
double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
double diff[NMAX];
double  epsilon = 0.0001;
int vector_n;
double ck, L_val, p0, p1;
double xstar[2], x0[2], value;

void subtract_vector(double dest[], double source[], int n)
{
	int i;
	for (i = 0; i < n; i++)
	     dest[i] = dest[i] - source[i];
} 

double falpha(double alpha)
{
	int i;
	printf("vector_n = %d, alpha = %lf\n", vector_n, alpha);
	for (i = 0; i < vector_n; i++)
		gtemp[i] = xnm1[i] - alpha * grad_vector[i];

	for (i = 0; i < vector_n; i++)
		printf("gtemp[%d] = %lf,  xnm1[%d] = %lf, grad_vector[%d] = %lf\n", i, gtemp[i], i, xnm1[i], i, grad_vector[i]);

	return  objective_function(gtemp);
} 
double golden(double (*fp)(double), double x1, double x3, double eps)
{
	double x2, fx2, fx3, x4, fx4;
	double phi = 1.618;
	double phi1 = 2.618;

	x2 = x1 + (x3 - x1) / phi1;
	fx2 = (*fp)(x2);
	x4 = x1 + (x3 - x1) / phi;
	fx4 = (*fp)(x4);

	do {

		if (fx2 > fx4)
		{
			x1 = x2;
			x2 = x4;
			fx2 = fx4;
			x4 = x1 + (x3 - x1) / phi;
			fx4 = (*fp)(x4);
		}
		else
		{
			x3 = x4;
			x4 = x2;
			fx4 = fx2;
			x2 = x1 + (x3 - x1) / phi1;
			fx2 = (*fp)(x2);
		}
	} while ((x3 - x1) > eps);

	return ((x1 + x3) / 2);
}
int vector_convergence_test(double arr[], int n, double epsilon)
{
	int i;
	for (i = 0; i < n; i++)
		printf("arr[%d] = %lf\n", i, arr[i]);

	for (i = 0; i < n; i++)
		if (fabs(arr[i]) > epsilon)
			return 0;
	return 1;
} 
void copy_vector(double dest[], double source[], int n)
{
	int i;
	for (i = 0; i < n; i++)
		dest[i] = source[i];
} 
void find_initial_alphas(double (*falpha)(double), double* alpha_2)
{
	int going_down_flag;
	double falpha1, falpha2, alpha1, alpha2;

	falpha1 = (*falpha)(0.0);

	alpha1 = 0.0;
	alpha2 = 0.0009765625; 

	going_down_flag = 1;

	while (going_down_flag == 1)
	{
		falpha2 = (*falpha)(alpha2);
		if (falpha2 >= falpha1)
			going_down_flag = 0;
		else
		{
			alpha1 = alpha2;
			falpha1 = falpha2;
			alpha2 = 2.0 * alpha2;
		}  
	} 
	printf("alpha1 = %lf, alpha2 = %lf\n", alpha1, alpha2);
	*alpha_2 = alpha2;
} 
void steepest(double xn[], double x0[], int n, FUN_PTR f, GRAD_FUN_PTR grad, double epsilon, VECTOR_CONVERGENCE_TEST v)
{
	double temp, alpha_1, alpha_2, alpha_k;
	int i;
	vector_n = n;
	copy_vector(xn, x0, n);
	grad(grad_vector, x0);
	objective_function = f;
	copy_vector(xnm1, xn, n);
	grad(grad_vector, xnm1);
	copy_vector(diff, xn, n);

	while ((v(diff, n, 0.001) == 0) || (v(grad_vector, n, 0.001) == 0))
	{

		find_initial_alphas(falpha, &alpha_2);

		alpha_k = golden(falpha, 0.0, alpha_2, epsilon);

		printf("alpha_k = %lf\n", alpha_k);

		for (i = 0; i < n; i++)
			diff[i] = alpha_k * grad_vector[i];
		for (i = 0; i < n; i++)
			xn[i] = xnm1[i] - diff[i];

		printf("xn:\n");
		for (i = 0; i < n; i++)
			printf(" xn[%d] = %lf\n", i, xn[i]);

		copy_vector(xnm1, xn, n);
		grad(grad_vector, xn);
	} 
}  
double h0(double x[])
{
	return  (x[0] * x[0] +
		x[1] * x[1] - 1.0);
} 
double h00(double x[])
{
	return  (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1.0);
}

double h1(double x[])
{
	return  (x[0] * x[1] + x[0] * x[2] + x[1] * x[2] - 1.0);
}
double f(double x[])
{
	double h_val;
	h_val = h0(x) * h0(x);
	return ((x[0] + 2 * x[1]) + ck * h_val);
} 
double f2(double x[])
{
	double h0_val, h1_val;
	h0_val = h0(x) * h0(x);
	h1_val = h1(x) * h1(x);
	return ((x[0] + 2 * x[1] + 3 * x[2]) + ck * (h0_val + h1_val));
}
double calc(double x[])
{
	return ((x[0] + 2 * x[1]));
}
double calc2(double x[])
{
	return (x[0] + 2 * x[1] + 3 * x[2]);
}
double approx_partial_derivative(double (*calc)(double x[]),
	int i, double x[])
{
	double temp1, temp2, xi_orig, result, h;
	double eps_const = 1048576.0;

	xi_orig = x[i];
	h = x[i] / eps_const;

	x[i] = xi_orig + h;

	temp1 = (*calc)(x);

	x[i] = xi_orig - h;

	temp2 = (*calc)(x);

	result = (temp1 - temp2) / (2 * h);

	x[i] = xi_orig;

	return result;

} 
void approx_g(double grad[], double x[])
{
	int i, j;

	for (i = 0; i < vector_n; i++)
		grad[i] = approx_partial_derivative(f, i, x);

} 
void affine_c_n(void)
{
	if ((p1 - p0) < 0.005) {
		ck = ck * 3 / 4; 
		if (L_val < 0.005)
			ck = ck * 7 / 6; 
		steepest(xstar, x0, 2, f, approx_g, 0.0001, vector_convergence_test);
		L_val = sqrt((x0[0] - xstar[0]) * (x0[0] - xstar[0]) + (x0[1] - xstar[1]) * (x0[1] - xstar[1]));
	}

}
void main()
{

	///////////
	double xstar[3], x0[3], value, L_val;

	x0[0] = -0.5;
	x0[1] = -0.5;
	x0[2] = -0.5;
	ck = 1.0;
	steepest(xstar, x0, 3, f2, approx_g, 0.004, vector_convergence_test);
	copy_vector(x0, xstar, 3);

	do {
		ck *= 2.0;
		steepest(xstar, x0, 3, f2, approx_g, 0.004, vector_convergence_test);
		L_val = sqrt((x0[0] - xstar[0]) * (x0[0] - xstar[0]) + (x0[1] - xstar[1]) * (x0[1] - xstar[1]) + (x0[2] - xstar[2]) * (x0[2] - xstar[2]));
		copy_vector(x0, xstar, 3);

	} while (L_val >= 0.005);



	printf("\n\n  optimal solution:\n   xstar[0] = %lf,\n   xstar[1] = %lf\n,  xstar[2] = % lf\n", xstar[0], xstar[1], xstar[2]);
	printf("Ck = %1f2\n", ck);
	printf("\n\n   optimal value  = %lf\n", calc2(xstar));



	///////////

	printf("***************************************************************\n\n\n");
	///////////
	x0[0] = -0.5;
	x0[1] = -0.5;
	ck = 1.0;
	steepest(xstar, x0, 2, f, approx_g, 0.0001, vector_convergence_test);
	copy_vector(x0, xstar, 2);
	do {
		ck *= 2.0;
		steepest(xstar, x0, 2, f, approx_g, 0.0001, vector_convergence_test);
		L_val = sqrt((x0[0] - xstar[0]) * (x0[0] - xstar[0]) + (x0[1] - xstar[1]) * (x0[1] - xstar[1]));
		copy_vector(x0, xstar, 2);

	} while (L_val >= 0.0001);
	printf("\n\noptimal solution:\n   xstar[0] = %lf,\n   xstar[1] = %lf\n", xstar[0], xstar[1]);
	printf("\n\nIn degrees: xstar[0] = %lf\n,  xstar[1] = %lf\n",
		xstar[0] * 180.0 / M_PI, xstar[1] * 180.0 / M_PI);
	printf("Ck = %1f\n", ck);
	printf("\n\noptimal value  = %lf\n", calc(xstar));
	///////////
	printf("***************************************************************\n\n\n");
	x0[0] = -0.5;
	x0[1] = -0.5;
	ck = 1.0;
	steepest(xstar, x0, 2, f, approx_g, 0.0001, vector_convergence_test);
	copy_vector(x0, xstar, 2);
	p1 = f(x0);
	do {
		ck = ck * 2;
		steepest(xstar, x0, 2, f, approx_g, 0.0001, vector_convergence_test);
		L_val = sqrt((x0[0] - xstar[0]) * (x0[0] - xstar[0]) + (x0[1] - xstar[1]) * (x0[1] - xstar[1]));
		copy_vector(x0, xstar, 2);
		p0 = p1;
		p1 = f(xstar);
	} while (L_val >= 0.005);
	affine_c_n();   
	printf("\n\n optimal solution:\n   xstar[0] = %lf,\n   xstar[1] = %lf\n", xstar[0], xstar[1]);
	printf("\nCq = %1f", ck);
	printf("\n\n optimal value  = %lf\n", calc(xstar));
}