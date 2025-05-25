#include "mathalg.h"

Matrix3ds Log(const Matrix3ds& p, const Matrix3ds& X)
{
	double l[3], s[3];
	Vector3d u[3], v[3];

	// evaluate eigen-decomposition of p
	p.eigen(l, u);
	Matrix3d U(u[0], u[1], u[2]);
	Matrix3dd L(l[0], l[1], l[2]);
	Matrix3dd rootL(sqrt(l[0]), sqrt(l[1]), sqrt(l[2]));

	Matrix3d G = U * rootL;
	Matrix3d Gi = G.inverse();

	Matrix3ds Y = (Gi * X*Gi.transpose()).sym();

	Y.eigen(s, v);
	Matrix3d V = Matrix3d(v[0], v[1], v[2]);

	Matrix3d GV = G * V;

	Matrix3dd logS(log(s[0]), log(s[1]), log(s[2]));

	Matrix3d LogX = (GV)*logS*(GV.transpose());

	return LogX.sym();
}

Matrix3ds Exp(const Matrix3ds& p, const Matrix3ds& X)
{
	double l[3], s[3];
	Vector3d u[3], v[3];

	// evaluate eigen-decomposition of p
	p.eigen(l, u);
	Matrix3d U(u[0], u[1], u[2]);
	Matrix3dd L(l[0], l[1], l[2]);
	Matrix3dd rootL(sqrt(l[0]), sqrt(l[1]), sqrt(l[2]));

	Matrix3d G = U * rootL;
	Matrix3d Gi = G.inverse();

	Matrix3ds Y = (Gi * X*Gi.transpose()).sym();

	Y.eigen(s, v);
	Matrix3d V = Matrix3d(v[0], v[1], v[2]);

	Matrix3d GV = G * V;

	Matrix3dd expS(exp(s[0]), exp(s[1]), exp(s[2]));

	Matrix3d ExpX = (GV)*expS*(GV.transpose());

	return ExpX.sym();
}

Matrix3ds weightedAverageStructureTensor(Matrix3ds* d, double* w, int n)
{
	const double eps = 1.0e-9;

	Matrix3ds mu(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
	double tau = 1.0;

	double normXi = 0.0, normXip = 0.0;
	Matrix3ds Xi, Xip;
	int nc = 0;
	do
	{
		Xip = Xi;
		normXip = normXi;

		Xi = weightedAverage<Matrix3ds>(d, w, n, [&](const Matrix3ds& a) {
			return Log(mu, a); 
		});

		mu = Exp(mu*tau, Xi);

		normXi = Xi.norm();

		if ((nc != 0) && (normXi > normXip))
		{
			Xi = Xip;
			tau *= 0.5;
		}
		nc++;
	} 
	while (normXi > eps);

	return mu;
}
