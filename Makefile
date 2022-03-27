# Makefile
#     Copyright (C) 2022  Rahul Dhodapkar
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

docs:
	echo "roxygen2::roxygenise()" | R --no-save
test:
	echo "devtools::test()" | R --no-save
build:
	echo "devtools::load_all()" | R --no-save
check:
	echo "devtools::check()" | R --no-save
release:
	echo "devtools::release()" | R --no-save
install:
	echo "devtools::install('../puddlr')" | R --no-save
