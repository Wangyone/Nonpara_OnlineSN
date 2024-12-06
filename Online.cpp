#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec statMarginal(arma::mat x, int T ,double alphaV,arma::vec vahaY,String kern,double delta) {  //x 时间序列，T 历史数据
	int nY = x.n_rows;    // N总样本量
	int d  = x.n_cols;     // 维数
	arma::mat Y=x;
	arma::mat Wjk(nY, nY);			  //距离矩阵Wjk
	
	int Tm = T ; 				//——> T-m+1
	arma::vec tstat(nY - Tm);  // 统计量的个数为N-T个
	
	if(kern=="energy"){
		for (int j = 0; j < nY; j++) {
		for (int k = 0; k < nY; k++) {
		Wjk(j,k) = pow( sum(pow(Y.row(j)-Y.row(k),2)) , alphaV/2); 
			}
			}   //距离矩阵--对应的平方和||Yj-Yk||^(alphaV/2)
		int Tmt = 0;              
		
	for (int tt = 0; tt < (nY - Tm); tt++) {  //搜寻变点
			Tmt = Tm + tt + 1;  //cpt=T-m+(1+tt)——> T+t-m
			//submat(0,0,Tm-1,Tm-1)子矩阵的起始行索引为0，起始列索引为0
			//结束行索引为Tm-1，结束列索引为Tm-1
			tstat(tt) = 2*sum(sum(Wjk.submat(0,0,Tmt-1,Tm-1)))/Tm/Tmt- 
				sum(sum(Wjk.submat(0,0,Tm-1,Tm-1)))/pow(Tm,2)-  
				sum(sum(Wjk.submat(0,0,Tmt-1,Tmt-1)))/pow(Tmt,2);
			tstat(tt) =vahaY(tt)*tstat(tt);
			}
			}
	else if(kern == "gauss"){  
		for (int j = 0; j < nY; j++) {
		for (int k = 0; k < nY; k++) {
		Wjk(j,k) = exp( -sum(pow(Y.row(j)-Y.row(k),2))/(4*delta)); 
			}
			}   //距离矩阵
		
		int Tmt = 0;              //cpt  nY=N-m+1
		
		for (int tt = 0; tt < (nY - Tm); tt++) {  //搜寻变点
			Tmt = Tm + tt + 1;  //cpt=T-m+1+(1+tt)——> T+t-m+1
			
			//submat(0,0,Tm-1,Tm-1)子矩阵的起始行索引为0，起始列索引为0
			//结束行索引为Tm-1，结束列索引为Tm-1
			tstat(tt) = -2*sum(sum(Wjk.submat(0,0,Tmt-1,Tm-1)))/Tm/Tmt+
				sum(sum(Wjk.submat(0,0,Tm-1,Tm-1)))/pow(Tm,2)+  
				sum(sum(Wjk.submat(0,0,Tmt-1,Tmt-1)))/pow(Tmt,2);
			
			tstat(tt) =vahaY(tt)*tstat(tt)*pow(datum::pi/delta, d); 
		}	   //使用arma::datum::pi来获取圆周率的值
	} 
	else if(kern == "quadexp"){  
		//推荐的核函数
		double d1=2*delta;
		double d2=2*d1;
		for (int j = 0; j < nY; j++) {
			for (int k = 0; k < nY; k++) {
				arma::vec v2=-pow(Y.row(j)-Y.row(k),2).t();
				arma::vec v3=exp(v2/d2)%(d1+v2)/d1;  // % schur product对应元素相乘
				// 将向量 v3 中的所有元素相乘
				Wjk(j,k) =prod(v3);
			}
		}   //距离矩阵 
		
		int Tmt = 0;              
		for (int tt = 0; tt < (nY - Tm); tt++) {  //搜寻变点
			Tmt = Tm + tt + 1;  //cpt=T-m+1+(1+tt)——> T+t-m+1
			
			//submat(0,0,Tm-1,Tm-1)子矩阵的起始行索引为0，起始列索引为0
			//结束行索引为Tm-1，结束列索引为Tm-1
			tstat(tt) = 2*sum(sum(Wjk.submat(0,0,Tmt-1,Tm-1)))/Tm/Tmt-
				sum(sum(Wjk.submat(0,0,Tm-1,Tm-1)))/pow(Tm,2)-  
				sum(sum(Wjk.submat(0,0,Tmt-1,Tmt-1)))/pow(Tmt,2);
			
			tstat(tt) =-1*vahaY(tt)*tstat(tt); 
		}	 
	}
	return tstat;
}


