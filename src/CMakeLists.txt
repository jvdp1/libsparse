# Based on test-drive.
# SPDX-Identifier: MIT
#
# Licensed under either of MIT license; you may not use this file except
# in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/dgtrsm.f"
  "${dir}/modcommon.f90"
  "${dir}/modhash.f90"
  "${dir}/modrandom.f90"
  "${dir}/modspainv.f90"
  "${dir}/modsparse_coo.f90"
  "${dir}/modsparse_crs64.f90"
  "${dir}/modsparse_crs.f90"
  "${dir}/modsparse.f90"
  "${dir}/modsparse_gen.f90"
  "${dir}/modsparse_ll.f90"
  "${dir}/modsparse_metisgraph.f90"
  "${dir}/modsparse_mkl.f90"
  "${dir}/modvariablepardiso.f90"
  "${dir}/sgtrsm.f"
  "${dir}/smbfct.f"
)

if(CMAKE_METIS_LIB)
  list(
    APPEND srcs
    "${dir}/modmetis.f90"
  )
endif()

set(srcs "${srcs}" PARENT_SCOPE)
