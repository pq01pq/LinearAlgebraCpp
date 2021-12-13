#include "linalg_allocate.h"

namespace linalg {
	Allocatorr::Allocatorr(Allocatablee& target, const size_t sequence)
		: target(target), sequence(sequence)
	{
	}
	Allocatorr::~Allocatorr()
	{
	}
	Allocatorr& Allocatorr::operator,(const double value)
	{
		target.allocate(sequence, value);
		sequence++;
		return *this;
	}
}