/*********************************************************************
 * \file   FEAppConfig.cpp
 * \brief
 *
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#include "app/FEAppConfig.h"

FEAppConfig::FEAppConfig()
{
    m_noutput = 1;
    m_bRunFile = true;
    Defaults();
}

void FEAppConfig::SetOutputLevel(int n)
{
    m_noutput = n;
}

bool FEAppConfig::ReadConfig(std::string filepath)
{
    return true;
}

void FEAppConfig::Defaults()
{
    m_printParams = -1;
    m_bshowErrors = true;
}
