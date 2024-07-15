#include "fieldmanager.h"

namespace fem
{
	FieldManager :: ~FieldManager(){ }

	void
		FieldManager::registerField(std::shared_ptr< Field >eField, FieldType key)
	{
		if (this->externalFields.find(key) == this->externalFields.end()) {
			this->externalFields.insert({ key, eField });
		}
		else {
			this->externalFields.at(key) = eField;
		}
		/// could be replaced in future by 
		/// this->externalFields.insert_or_assign({key, eField});
	}




	bool
		FieldManager::isFieldRegistered(FieldType key)
	{
		return (this->externalFields.find(key) != this->externalFields.end());
	}

	std::shared_ptr< Field >
		FieldManager::giveField(FieldType key)
	{
		auto i = this->externalFields.find(key);
		if (i == this->externalFields.end()) {
			std::shared_ptr< Field >p; // std::shared_ptr<Field> p(nullptr);
			return p;
		}

		return i->second;
	}

	void
		FieldManager::unregisterField(FieldType key)
	{
		auto i = this->externalFields.find(key);
		if (i == this->externalFields.end()) {
			return;
		}

		this->externalFields.erase(i);
	}

	std::vector<FieldType>
		FieldManager::giveRegisteredKeys()
	{
		std::vector<FieldType> ret;
		for (const auto& keyField : this->externalFields) ret.push_back(keyField.first);
		return ret;
	}

} // end namespace oofem
