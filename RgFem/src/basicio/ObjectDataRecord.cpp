#include "ObjectDataRecord.h"

//-----------------------------------------------------------------------------
ObjectDataRecord::ObjectDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_RB) 
{

}

//-----------------------------------------------------------------------------
void ObjectDataRecord::SetData(const char* szexpr)
{
	/*char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FELogObjectData* pdata = fecore_new<FELogObjectData>(sz, GetFEModel());
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);*/
}

//-----------------------------------------------------------------------------
double ObjectDataRecord::Evaluate(int item, int ndata)
{
	//FEMechModel* fem = dynamic_cast<FEMechModel*>(GetFEModel());

	//FEMesh& mesh = fem->GetMesh();
	//int nrb = item - 1;
	//if ((nrb < 0) || (nrb >= fem->Materials())) return 0;

	//double val = 0;

	//// find the rigid body that has this material
	//int NRB = fem->RigidBodies();
	//for (int i=0; i<NRB; ++i)
	//{
	//	FERigidBody& obj = *fem->GetRigidBody(i);
	//	if (obj.GetMaterialID() == nrb) return m_Data[ndata]->value(obj);
	//}

	return 0;
}

//-----------------------------------------------------------------------------
int ObjectDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
void ObjectDataRecord::SelectAllItems()
{
    /*FEMechModel* fem = dynamic_cast<FEMechModel*>(GetFEModel());

    int n = 0, i;
    for (i=0; i<fem->Materials(); ++i)
    {
        FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem->GetMaterial(i));
        if (pm) ++n;
    }

    if (n > 0)
    {
        m_item.resize(n);
        n = 0;
        for (i=0; i<fem->Materials(); ++i)
        {
            FERigidMaterial* pm  = dynamic_cast<FERigidMaterial*>(fem->GetMaterial(i));
            if (pm)
            {
                m_item[n++] = i+1;
            }
        }
    }*/
}
