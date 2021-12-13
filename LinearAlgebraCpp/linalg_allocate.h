#pragma once

namespace linalg {

	class Allocatorr;
	class Allocatablee;

	/*
	* Allocator class is used for allocating values into linear algebra containers.
	* It allocates values when it encounters operator '<<' and ',' with increasing sequence.
	*/
	class Allocatorr {
		friend class Allocatablee;
	public:
		Allocatorr(Allocatablee& target, const size_t sequence);
		~Allocatorr();
		Allocatorr& operator,(const double value);

	private:
		Allocatablee& target;
		size_t sequence;
	};

	// Mix-in class of linear algebra containers(== Java interface)
	class Allocatablee {
		friend class Allocatorr;
	protected:
		virtual void allocate(const size_t sequence, const double value) = 0;
	};
}