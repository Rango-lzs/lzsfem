/*****************************************************************//**
 * \file   Material.h
 * \brief  
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef MATERIAL_H
#define MATERIAL_H

#include<iostream>
#include<string>
#include"../algebra/Matrix.h"
#include"../algebra/Vector.h"
#include"../tensors/Tensor.h"

/*
* Task: 
*	give the matrial constitutive matrix which related the the stress and strain status
*	ʵ�ֲ��ϱ���ģ�ͣ����ϸնȣ�Ӧ������
*   ��ʷ�����Ĵ洢
*/

class Material
{
	public:
		Material()  = default;

	    virtual Matrix<double>& constitutiveMatrix() = 0;
	
		virtual std::string matType() = 0;
			
		virtual void assembleTensors(Vector<double>& v, Tensor& strains, Tensor& stresses)=0;

};


#endif