#ifndef PROGRESS_INDICATOR_STRATEGY_HPP
#define PROGRESS_INDICATOR_STRATEGY_HPP

#include <string>

/*
This class implements a Strategy pattern to enable the possibility of tweaking
the progress indicator code without messing with the analysis code These
routines are then used to update the GUI
*/
class ProgressIndicatorStrategy {
protected:
	std::string m_current_section_name;
	size_t m_progress_limit;
	size_t m_current_progress;

public:
	ProgressIndicatorStrategy() : m_progress_limit(0), m_current_progress(0) {}
	/**
							Marks the begining of a new progress section
							@param	section name
							**/
	virtual void markSectionStart(std::string) {};

	/**
							Sets the iterations range that the current section
	   must go through
							**/
	virtual void markSectionLimit(size_t) {};


	/**
							Increments the current iterator
							**/
	virtual void markSectionIterationIncrement() {};

	/**
							Marks the end of the current progress section
							**/
	virtual void markSectionEnd() {};

	/**
							Sets the current progress
							@param	progress
							**/
	virtual void markProgress(size_t) {};

	virtual void message(std::string) {};

	virtual void error(std::string) {};

	/**
							Announces the end of the entire process
							**/
	virtual void markFinish() {};
};

#endif
