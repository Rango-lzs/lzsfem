/*********************************************************************
 * \file   FEAppConfig.h
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once
#include "femcore/fem_export.h"

class FEM_EXPORT FEAppConfig
{
public:
    FEAppConfig();
    void Defaults();
    void SetOutputLevel(int n);
    bool ReadConfig(std::string filepath);

public:
    int m_printParams;
    int m_noutput;
    bool m_bshowErrors;
};
