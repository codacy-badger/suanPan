﻿/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifndef SPARSE_FORM_HPP
#define SPARSE_FORM_HPP

#include <suanPan.h>
#include <armadillo>
#include <vector>

#ifdef SUANPAN_MT
#include <tbb/parallel_sort.h>
#define suanpan_sort tbb::parallel_sort
#else
#define suanpan_sort std::sort
#endif

using namespace arma;

#ifdef SUANPAN_MAGMA
using index_t = magma_index_t;
#else
using index_t = uword;
#endif

template<typename T, typename FT> class sparse_form {
protected:
	const T bin = 0.; // bin for out of bound elements

	virtual void copy_memory(index_t, const index_t*, const index_t*, const T*) = 0;
public:
	typedef FT form_type;

	const index_t n_rows = 0; // number of rows
	const index_t n_cols = 0; // number of cols
	const index_t n_elem = 0; // maximum number of elements
	const index_t c_size = 0; // current number of valid elements

	sparse_form() = default;
	sparse_form(index_t, index_t, index_t = 0);

	virtual ~sparse_form() = default;

	sparse_form(const sparse_form&) = delete;            // copy ctor
	sparse_form(sparse_form&&) = delete;                 // move ctor
	sparse_form& operator=(const sparse_form&) = delete; // copy assignment
	sparse_form& operator=(sparse_form&&) = delete;      // move assignment

	virtual bool is_empty() const;

	virtual void reset() const = 0;
	virtual void zeros() const = 0;

	virtual T max() const = 0;

	virtual bool init() = 0;
	virtual bool init(index_t) = 0;
	virtual bool init(index_t, index_t, index_t) = 0;
	virtual bool resize() = 0;
	virtual bool resize(index_t) = 0;
	virtual bool resize(index_t, index_t, index_t) = 0;

	virtual void print() const;
	virtual void spy();
};

template<typename T, typename FT> sparse_form<T, FT>::sparse_form(const index_t in_rows, const index_t in_cols, const index_t in_elem)
	: n_rows(in_rows)
	, n_cols(in_cols)
	, n_elem(in_elem) {}

template<typename T, typename FT> bool sparse_form<T, FT>::is_empty() const { return c_size == 0; }

template<typename T, typename FT> void sparse_form<T, FT>::print() const {}

template<typename T, typename FT> void sparse_form<T, FT>::spy() { throw invalid_argument("not supported"); }

#endif
