/*
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

#ifndef BAMFILE_HPP
#define BAMFILE_HPP

#include <string>
#include "GenomicRegion.hpp"
#include <htslib/sam.h>

class BAMFile
{
public:
    BAMFile();
    BAMFile(const std::string &filename,
            const std::string &mode = "");
    ~BAMFile();

    void
    open(const std::string &filename,
         const std::string &mode = "");

    void
    close();

    bool
    good() const {return GOOD;}

    friend BAMFile&
    operator>>(BAMFile& in, SimpleGenomicRegion& region);

    friend BAMFile&
    operator>>(BAMFile& in, GenomicRegion& region);

private:
    std::string filename;
    std::string mode;

    htsFile*  file_handler;
    bam_hdr_t *hdr;
    bool GOOD;
};

#endif

