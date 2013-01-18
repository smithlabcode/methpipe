d/*
 * Copyright (C) 2012 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "smithlab_utils.hpp"

using std::vector;
using std::string;
using std::endl;
using std::cerr;
using std::getline;

static void
convert_to_stringstream(const string &infile, std::stringstream &ss)
{
    std::ifstream in(infile.c_str());
    
    while (!in.eof())
    {
        string str;
        getline(in, str);
        const size_t comment_start = str.find("#");
        if (comment_start != string::npos)
            str.erase(comment_start);
        str = smithlab::strip(str);
        if (!str.empty())
            ss << str << endl;
    }
    in.close();
}

template<class Distro_Type> void
read_param_file(const std::string &infile, size_t &n,
                std::vector<std::vector<double> > &trans,
                std::vector<Distro_Type> &emissions,
                std::vector<Distro_Type> &durations);
{
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    convert_to_stringstream(infile, ss);
    
    ss >> n; 
    string tmp_str;
    getline(ss, tmp_str);

    emissions.clear();
    for (size_t i = 0; i < n; ++i)
    {
        string tmp_str;
        getline(ss, tmp_str);
        emissions.push_back(Distro_Type(tmp_str));
    }

    durations.clear();
    for (size_t i = 0; i < n; ++i)
    {
        string tmp_str;
        getline(ss, tmp_str);
        durations.push_back(Distro_Type(tmp_str));
    }
    
    trans.resize(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            ss >> trans[i][j];
}

template <class Distro_Type> void
write_param_file(const std::string &outfile,  const size_t &n,
                 const std::vector<std::vector<double> > &trans,
                 const std::vector<Distro_Type> &emissions,
                 const std::vector<Distro_Type> &durations);
{
    std::ofstream out(outfile.c_str());
    
    out << "# number of states" << endl;
    out << n << endl;
    
    out << "\n# emmission distributions" << endl;
    std::copy(emissions.begin(), emissions.end(),
              std::ostream_iterator<Distro_Type>(out, "\n"));

    out << "\n# duration distributions" << endl;
    std::copy(durations.begin(), durations.end(),
              std::ostream_iterator<Distro_Type>(out, "\n"));
    
    out << "\n# state transition probabilities" << endl;
    for (size_t i = 0; i < n; ++i)
    {
        copy(trans[i].begin(), trans[i].end(),
             std::ostream_iterator<double>(out, "\t"));
        out << endl;
    }
}

