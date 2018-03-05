// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This is the unit test for the black oil fluid system
 *
 * This test requires the presence of opm-parser.
 */
#include "config.h"

#if !HAVE_ECL_INPUT
#error "The test for the black oil fluid system classes requires ecl input support in opm-common"
#endif

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <type_traits>
#include <cmath>

// values of strings based on the SPE1 and NORNE cases of opm-data.
static const char* deckString1 =
    "RUNSPEC\n"
    "\n"
    "DIMENS\n"
    "   10 10 3 /\n"
    "\n"
    "METRIC\n"
    "\n"
    "TABDIMS\n"
    " * 2 /\n"
    "\n"
    "OIL\n"
    "GAS\n"
    "WATER\n"
    "\n"
    "DISGAS\n"
    "VAPOIL\n"
    "\n"
    "GRID\n"
    "\n"
    "DX\n"
    "   	300*1000 /\n"
    "DY\n"
    "	300*1000 /\n"
    "DZ\n"
    "	100*20 100*30 100*50 /\n"
    "\n"
    "TOPS\n"
    "	100*1234 /\n"
    "\n"
    "PROPS\n"
    "\n"
    "-- Copyright (C) 2015 Statoil\n"
    "\n"
    "PVTG\n"
    "\n"
    "-- PRESSURE       RSG        B-GAS     VISCOSITY\n"
    "--   BAR                                 (CP)\n"
    "\n"
    "     50.00    0.00000497   0.024958     0.01441\n"
    "              0.00000248   0.024958     0.01440\n"
    "              0.00000000   0.024958     0.01440 /\n"
    "\n"
    "     70.00    0.00000521   0.017639     0.01491\n"
    "              0.00000261   0.017641     0.01490\n"
    "              0.00000000   0.017643     0.01490 /\n"
    "\n"
    "     90.00    0.00000627   0.013608     0.01547\n"
    "              0.00000313   0.013611     0.01546\n"
    "              0.00000000   0.013615     0.01544 /\n"
    "\n"
    "    110.00    0.00000798   0.011072     0.01609\n"
    "              0.00000399   0.011076     0.01607\n"
    "              0.00000000   0.011081     0.01605 /\n"
    "\n"
    "    130.00    0.00001041   0.009340     0.01677\n"
    "              0.00000520   0.009346     0.01674\n"
    "              0.00000000   0.009352     0.01671 /\n"
    "\n"
    "    150.00    0.00001365   0.008092     0.01752\n"
    "              0.00000683   0.008099     0.01748\n"
    "              0.00000000   0.008106     0.01743 /\n"
    "\n"
    "    170.00    0.00001786   0.007156     0.01834\n"
    "              0.00000893   0.007164     0.01827\n"
    "              0.00000000   0.007172     0.01819 /\n"
    "\n"
    "    190.00    0.00002316   0.006433     0.01923\n"
    "              0.00001158   0.006442     0.01912\n"
    "              0.00000000   0.006451     0.01900 /\n"
    "\n"
    "    210.00    0.00002972   0.005861     0.02019\n"
    "              0.00001486   0.005871     0.02001\n"
    "              0.00000000   0.005881     0.01984 /\n"
    "\n"
    "    230.00    0.00003767   0.005402     0.02121\n"
    "              0.00001883   0.005412     0.02095\n"
    "              0.00000000   0.005422     0.02071 /\n"
    "\n"
    "    250.80    0.00004756   0.005013     0.02234\n"
    "              0.00002378   0.005022     0.02197\n"
    "              0.00000000   0.005032     0.02162 /\n"
    "\n"
    "    268.42    0.00005757   0.004737     0.02335\n"
    "              0.00002878   0.004746     0.02287\n"
    "              0.00000000   0.004754     0.02240 /\n"
    "\n"
    "    285.33    0.00006853   0.004511     0.02438\n"
    "              0.00003427   0.004518     0.02375\n"
    "              0.00000000   0.004525     0.02315 /\n"
    "\n"
    "    301.59    0.00008041   0.004323     0.02542\n"
    "              0.00004020   0.004327     0.02463\n"
    "              0.00000000   0.004332     0.02387 /\n"
    "\n"
    "    317.23    0.00009313   0.004165     0.02648\n"
    "              0.00004657   0.004166     0.02549\n"
    "              0.00000000   0.004169     0.02456 /\n"
    "\n"
    "    332.29    0.00010668   0.004031     0.02755\n"
    "              0.00005334   0.004029     0.02634\n"
    "              0.00000000   0.004028     0.02522 /\n"
    "\n"
    "    346.80    0.00012100   0.003917     0.02863\n"
    "              0.00006050   0.003911     0.02719\n"
    "              0.00000000   0.003906     0.02585 /\n"
    "\n"
    "    360.80    0.00013607   0.003819     0.02974\n"
    "              0.00006803   0.003808     0.02803\n"
    "              0.00000000   0.003799     0.02645 /\n"
    "\n"
    "    374.31    0.00015188   0.003735     0.03087\n"
    "              0.00007594   0.003718     0.02887\n"
    "              0.00000000   0.003705     0.02703 /\n"
    "\n"
    "    387.36    0.00016843   0.003662     0.03202\n"
    "              0.00008421   0.003639     0.02970\n"
    "              0.00000000   0.003621     0.02758 /\n"
    "\n"
    "    399.99    0.00018571   0.003598     0.03320\n"
    "              0.00009286   0.003570     0.03053\n"
    "              0.00000000   0.003545     0.02810 /\n"
    "\n"
    "    412.21    0.00020375   0.003543     0.03442\n"
    "              0.00010188   0.003508     0.03137\n"
    "              0.00000000   0.003477     0.02861 /\n"
    "\n"
    "    424.05    0.00022256   0.003496     0.03566\n"
    "              0.00011128   0.003453     0.03220\n"
    "              0.00000000   0.003416     0.02909 /\n"
    "\n"
    "    435.53    0.00024218   0.003454     0.03695\n"
    "              0.00012109   0.003404     0.03305\n"
    "              0.00000000   0.003360     0.02956 /\n"
    "\n"
    "    446.68    0.00026266   0.003419     0.03828\n"
    "              0.00013133   0.003360     0.03390\n"
    "              0.00000000   0.003309     0.03000 /\n"
    "\n"
    "    457.51    0.00028404   0.003388     0.03967\n"
    "              0.00014202   0.003320     0.03477\n"
    "              0.00000000   0.003262     0.03043 /\n"
    "\n"
    "    468.04    0.00030639   0.003362     0.04110\n"
    "              0.00015319   0.003285     0.03566\n"
    "              0.00000000   0.003218     0.03085 /\n"
    "\n"
    "    478.30    0.00032980   0.003341     0.04261\n"
    "              0.00016490   0.003253     0.03656\n"
    "              0.00000000   0.003178     0.03125 /\n"
    "\n"
    "    488.30    0.00035436   0.003323     0.04418\n"
    "              0.00017718   0.003225     0.03749\n"
    "              0.00000000   0.003140     0.03164 /\n"
    "\n"
    "    498.06    0.00038020   0.003310     0.04583\n"
    "              0.00019010   0.003200     0.03845\n"
    "              0.00000000   0.003105     0.03202 /\n"
    "\n"
    "    507.59    0.00040745   0.003300     0.04758\n"
    "              0.00020373   0.003178     0.03944\n"
    "              0.00000000   0.003073     0.03238 /\n"
    "\n"
    "    516.92    0.00043630   0.003293     0.04943\n"
    "              0.00021815   0.003158     0.04048\n"
    "              0.00000000   0.003042     0.03273 /\n"
    "\n"
    "    526.06    0.00046694   0.003290     0.05141\n"
    "              0.00023347   0.003141     0.04156\n"
    "              0.00000000   0.003013     0.03308 /\n"
    "\n"
    "    535.02    0.00049963   0.003291     0.05353\n"
    "              0.00024981   0.003126     0.04271\n"
    "              0.00000000   0.002986     0.03342 /\n"
    "\n"
    "    543.83    0.00053469   0.003295     0.05582\n"
    "              0.00026734   0.003114     0.04393\n"
    "              0.00000000   0.002960     0.03374 /\n"
    "\n"
    "    552.49    0.00057251   0.003303     0.05830\n"
    "              0.00028625   0.003105     0.04523\n"
    "              0.00000000   0.002935     0.03407 /\n"
    "\n"
    "    561.04    0.00061359   0.003315     0.06103\n"
    "              0.00030679   0.003098     0.04664\n"
    "              0.00000000   0.002912     0.03438 /\n"
    "\n"
    "    569.48    0.00065855   0.003332     0.06405\n"
    "              0.00032928   0.003093     0.04818\n"
    "              0.00000000   0.002890     0.03469 /\n"
    "\n"
    "    577.82    0.00070820   0.003354     0.06744\n"
    "              0.00035410   0.003092     0.04988\n"
    "              0.00000000   0.002868     0.03500 /\n"
    "\n"
    "    586.09    0.00076355   0.003382     0.07127\n"
    "              0.00038178   0.003094     0.05178\n"
    "              0.00000000   0.002847     0.03530 /\n"
    "\n"
    "    594.29    0.00082592   0.003418     0.07567\n"
    "              0.00041296   0.003099     0.05394\n"
    "              0.00000000   0.002828     0.03560 /\n"
    "/\n"
    "\n"
    "-- PVT region 2 --\n"
    "\n"
    "     80.00    0.00000485   0.015151     0.01506\n"
    "              0.00000242   0.015154     0.01505\n"
    "              0.00000000   0.015157     0.01504 /\n"
    "\n"
    "\n"
    "    100.00    0.00000621   0.012032     0.01566\n"
    "              0.00000310   0.012036     0.01564\n"
    "              0.00000000   0.012040     0.01563 /\n"
    "\n"
    "\n"
    "    120.00    0.00000821   0.009980     0.01632\n"
    "              0.00000411   0.009985     0.01630\n"
    "              0.00000000   0.009990     0.01628 /\n"
    "\n"
    "\n"
    "    140.00    0.00001096   0.008537     0.01706\n"
    "              0.00000548   0.008544     0.01702\n"
    "              0.00000000   0.008550     0.01698 /\n"
    "\n"
    "\n"
    "    160.00    0.00001457   0.007476     0.01786\n"
    "              0.00000728   0.007484     0.01780\n"
    "              0.00000000   0.007492     0.01774 /\n"
    "\n"
    "\n"
    "    180.00    0.00001918   0.006669     0.01873\n"
    "              0.00000959   0.006678     0.01864\n"
    "              0.00000000   0.006687     0.01854 /\n"
    "\n"
    "\n"
    "    200.00    0.00002493   0.006038     0.01967\n"
    "              0.00001247   0.006048     0.01952\n"
    "              0.00000000   0.006058     0.01939 /\n"
    "\n"
    "\n"
    "    216.50    0.00003061   0.005616     0.02049\n"
    "              0.00001530   0.005626     0.02029\n"
    "              0.00000000   0.005636     0.02010 /\n"
    "/\n"
    "\n"
    "PVTO\n"
    "\n"
    "--  RSO    PRESSURE    B-OIL     VISCOSITY\n"
    "--          (BAR)                  (CP)\n"
    "\n"
    "   20.59     50.00    1.10615     1.180\n"
    "             75.00    1.10164     1.247\n"
    "            100.00    1.09744     1.315\n"
    "            125.00    1.09351     1.384\n"
    "            150.00    1.08984     1.453  /\n"
    "\n"
    "   28.19     70.00    1.12522     1.066\n"
    "             95.00    1.12047     1.124\n"
    "            120.00    1.11604     1.182\n"
    "            145.00    1.11191     1.241\n"
    "            170.00    1.10804     1.300  /\n"
    "\n"
    "   36.01     90.00    1.14458     0.964\n"
    "            115.00    1.13959     1.014\n"
    "            140.00    1.13494     1.064\n"
    "            165.00    1.13060     1.115\n"
    "            190.00    1.12653     1.166  /\n"
    "\n"
    "   44.09    110.00    1.16437     0.880\n"
    "            135.00    1.15915     0.924\n"
    "            160.00    1.15428     0.968\n"
    "            185.00    1.14973     1.012\n"
    "            210.00    1.14547     1.056  /\n"
    "\n"
    "   52.46    130.00    1.18467     0.805\n"
    "            155.00    1.17921     0.843\n"
    "            180.00    1.17413     0.882\n"
    "            205.00    1.16937     0.920\n"
    "            230.00    1.16491     0.959  /\n"
    "\n"
    "   61.13    150.00    1.20555     0.746\n"
    "            175.00    1.19985     0.780\n"
    "            200.00    1.19454     0.814\n"
    "            225.00    1.18958     0.849\n"
    "            250.00    1.18492     0.883  /\n"
    "\n"
    "   70.14    170.00    1.22704     0.698\n"
    "            195.00    1.22111     0.729\n"
    "            220.00    1.21558     0.759\n"
    "            245.00    1.21040     0.790\n"
    "            270.00    1.20555     0.821  /\n"
    "\n"
    "   79.50    190.00    1.24922     0.658\n"
    "            215.00    1.24305     0.686\n"
    "            240.00    1.23729     0.714\n"
    "            265.00    1.23190     0.742\n"
    "            290.00    1.22685     0.770  /\n"
    "\n"
    "   89.24    210.00    1.27214     0.637\n"
    "            235.00    1.26573     0.664\n"
    "            260.00    1.25974     0.693\n"
    "            285.00    1.25414     0.725\n"
    "            310.00    1.24888     0.760  /\n"
    "\n"
    "   99.39    230.00    1.29586     0.622\n"
    "            255.00    1.28921     0.641\n"
    "            280.00    1.28300     0.661\n"
    "            305.00    1.27718     0.680\n"
    "            330.00    1.27171     0.699  /\n"
    "\n"
    "  110.41    250.80    1.32148     0.610\n"
    "            275.80    1.31457     0.628\n"
    "            300.80    1.30812     0.647\n"
    "            325.80    1.30207     0.665\n"
    "            350.80    1.29638     0.682  /\n"
    "\n"
    "  120.32    268.42    1.34449     0.576\n"
    "            293.42    1.33735     0.593\n"
    "            318.42    1.33068     0.609\n"
    "            343.42    1.32442     0.626\n"
    "            368.42    1.31853     0.642  /\n"
    "\n"
    "  130.23    285.33    1.36737     0.5335\n"
    "            310.33    1.36001     0.5487\n"
    "            335.33    1.35313     0.5638\n"
    "            360.33    1.34667     0.5787\n"
    "            385.33    1.34059     0.5934  /\n"
    "\n"
    "  140.12    301.59    1.39015     0.4956\n"
    "            326.59    1.38257     0.5094\n"
    "            351.59    1.37548     0.5230\n"
    "            376.59    1.36882     0.5365\n"
    "            401.59    1.36255     0.5498  /\n"
    "\n"
    "  150.01    317.23    1.41282     0.4614\n"
    "            342.23    1.40503     0.4739\n"
    "            367.23    1.39773     0.4863\n"
    "            392.23    1.39088     0.4986\n"
    "            417.23    1.38443     0.5107  /\n"
    "\n"
    "  159.89    332.29    1.43539     0.43042\n"
    "            357.29    1.42739     0.44183\n"
    "            382.29    1.41990     0.45312\n"
    "            407.29    1.41286     0.46430\n"
    "            432.29    1.40622     0.47537 /\n"
    "\n"
    "  169.76    346.80    1.45788     0.41191\n"
    "            371.80    1.44967     0.42260\n"
    "            396.80    1.44198     0.43318\n"
    "            421.80    1.43475     0.44365\n"
    "            446.80    1.42794     0.45402 /\n"
    "\n"
    "  179.63    360.80    1.48028     0.39503\n"
    "            385.80    1.47187     0.40508\n"
    "            410.80    1.46398     0.41502\n"
    "            435.80    1.45657     0.42487\n"
    "            460.80    1.44958     0.43461 /\n"
    "\n"
    "  189.48    374.31    1.50260     0.37959\n"
    "            399.31    1.49399     0.38907\n"
    "            424.31    1.48591     0.39845\n"
    "            449.31    1.47832     0.40773\n"
    "            474.31    1.47116     0.41692 /\n"
    "\n"
    "  199.34    387.36    1.52484     0.36543\n"
    "            412.36    1.51603     0.37439\n"
    "            437.36    1.50777     0.38326\n"
    "            462.36    1.50000     0.39203\n"
    "            487.36    1.49267     0.40072 /\n"
    "\n"
    "  209.18    399.99    1.54700     0.35239\n"
    "            424.99    1.53800     0.36089\n"
    "            449.99    1.52956     0.36929\n"
    "            474.99    1.52161     0.37762\n"
    "            499.99    1.51411     0.38585 /\n"
    "\n"
    "  219.02    412.21    1.56910     0.34035\n"
    "            437.21    1.55991     0.34843\n"
    "            462.21    1.55128     0.35642\n"
    "            487.21    1.54316     0.36433\n"
    "            512.21    1.53549     0.37216 /\n"
    "\n"
    "  228.85    424.05    1.59112     0.32921\n"
    "            449.05    1.58174     0.33691\n"
    "            474.05    1.57294     0.34453\n"
    "            499.05    1.56464     0.35206\n"
    "            524.05    1.55681     0.35952 /\n"
    "\n"
    "  238.67    435.53    1.61307     0.31888\n"
    "            460.53    1.60351     0.32623\n"
    "            485.53    1.59453     0.33350\n"
    "            510.53    1.58606     0.34070\n"
    "            535.53    1.57807     0.34782 /\n"
    "\n"
    "  248.48    446.68    1.63496     0.30927\n"
    "            471.68    1.62522     0.31630\n"
    "            496.68    1.61606     0.32326\n"
    "            521.68    1.60743     0.33014\n"
    "            546.68    1.59927     0.33695 /\n"
    "\n"
    "  258.29    457.51    1.65678     0.30032\n"
    "            482.51    1.64686     0.30706\n"
    "            507.51    1.63753     0.31373\n"
    "            532.51    1.62873     0.32032\n"
    "            557.51    1.62042     0.32685 /\n"
    "\n"
    "  268.09    468.04    1.67853     0.29196\n"
    "            493.04    1.66843     0.29843\n"
    "            518.04    1.65893     0.30483\n"
    "            543.04    1.64997     0.31117\n"
    "            568.04    1.64150     0.31743 /\n"
    "\n"
    "  277.89    478.30    1.70022     0.28414\n"
    "            503.30    1.68994     0.29037\n"
    "            528.30    1.68028     0.29652\n"
    "            553.30    1.67116     0.30261\n"
    "            578.30    1.66253     0.30864 /\n"
    "\n"
    "  287.68    488.30    1.72184     0.27681\n"
    "            513.30    1.71139     0.28281\n"
    "            538.30    1.70156     0.28874\n"
    "            563.30    1.69228     0.29460\n"
    "            588.30    1.68350     0.30040 /\n"
    "\n"
    "  297.46    498.06    1.74339     0.26994\n"
    "            523.06    1.73277     0.27572\n"
    "            548.06    1.72278     0.28144\n"
    "            573.06    1.71334     0.28709\n"
    "            598.06    1.70442     0.29269 /\n"
    "\n"
    "  307.23    507.59    1.76487     0.26347\n"
    "            532.59    1.75409     0.26906\n"
    "            557.59    1.74393     0.27458\n"
    "            582.59    1.73434     0.28004\n"
    "            607.59    1.72527     0.28544 /\n"
    "\n"
    "  317.00    516.92    1.78628     0.25738\n"
    "            541.92    1.77533     0.26279\n"
    "            566.92    1.76502     0.26812\n"
    "            591.92    1.75528     0.27340\n"
    "            616.92    1.74606     0.27863 /\n"
    "\n"
    "  326.76    526.06    1.80761     0.25165\n"
    "            551.06    1.79651     0.25688\n"
    "            576.06    1.78604     0.26204\n"
    "            601.06    1.77615     0.26716\n"
    "            626.06    1.76679     0.27221 /\n"
    "\n"
    "  336.51    535.02    1.82887     0.24623\n"
    "            560.02    1.81761     0.25130\n"
    "            585.02    1.80699     0.25631\n"
    "            610.02    1.79696     0.26126\n"
    "            635.02    1.78746     0.26616 /\n"
    "\n"
    "  346.26    543.83    1.85005     0.24112\n"
    "            568.83    1.83864     0.24603\n"
    "            593.83    1.82787     0.25089\n"
    "            618.83    1.81770     0.25570\n"
    "            643.83    1.80806     0.26045 /\n"
    "\n"
    "  356.00    552.49    1.87115     0.23628\n"
    "            577.49    1.85959     0.24105\n"
    "            602.49    1.84868     0.24577\n"
    "            627.49    1.83836     0.25043\n"
    "            652.49    1.82858     0.25505 /\n"
    "\n"
    "  365.73    561.04    1.89217     0.23170\n"
    "            586.04    1.88046     0.23634\n"
    "            611.04    1.86940     0.24092\n"
    "            636.04    1.85895     0.24546\n"
    "            661.04    1.84904     0.24994 /\n"
    "\n"
    "  375.46    569.48    1.91309     0.22736\n"
    "            594.48    1.90124     0.23187\n"
    "            619.48    1.89004     0.23633\n"
    "            644.48    1.87946     0.24074\n"
    "            669.48    1.86942     0.24510 /\n"
    "\n"
    "  385.18    577.82    1.93391     0.22325\n"
    "            602.82    1.92192     0.22764\n"
    "            627.82    1.91060     0.23198\n"
    "            652.82    1.89988     0.23627\n"
    "            677.82    1.88971     0.24052 /\n"
    "\n"
    "  394.89    586.09    1.95464     0.21934\n"
    "            611.09    1.94252     0.22362\n"
    "            636.09    1.93106     0.22785\n"
    "            661.09    1.92021     0.23204\n"
    "            686.09    1.90993     0.23617 /\n"
    "\n"
    "  404.60    594.29    1.97527     0.21564\n"
    "            619.29    1.96301     0.21981\n"
    "            644.29    1.95143     0.22393\n"
    "            669.29    1.94046     0.22801\n"
    "            694.29    1.93005     0.23204 /\n"
    "/\n"
    "\n"
    "-- PVT region 2\n"
    "\n"
    "   32.91     80.00    1.13304     1.04537\n"
    "            114.00    1.12837     1.10009\n"
    "            148.00    1.12401     1.15521\n"
    "            182.00    1.11994     1.21063 /\n"
    "\n"
    "\n"
    "   40.99    100.00    1.15276     0.97219\n"
    "            134.00    1.14786     1.02086\n"
    "            168.00    1.14328     1.06983\n"
    "            202.00    1.13900     1.11901 /\n"
    "\n"
    "\n"
    "   49.36    120.00    1.17297     0.91124\n"
    "            154.00    1.16783     0.95496\n"
    "            188.00    1.16303     0.99891\n"
    "            222.00    1.15854     1.04301 /\n"
    "\n"
    "\n"
    "   58.04    140.00    1.19374     0.85942\n"
    "            174.00    1.18836     0.89902\n"
    "            208.00    1.18334     0.93878\n"
    "            242.00    1.17864     0.97864 /\n"
    "\n"
    "\n"
    "   67.04    160.00    1.21512     0.81456\n"
    "            194.00    1.20951     0.85065\n"
    "            228.00    1.20426     0.88686\n"
    "            262.00    1.19935     0.92313 /\n"
    "\n"
    "\n"
    "   76.39    180.00    1.23718     0.77508\n"
    "            214.00    1.23132     0.80815\n"
    "            248.00    1.22585     0.84130\n"
    "            282.00    1.22073     0.87448 /\n"
    "\n"
    "\n"
    "   86.11    200.00    1.25996     0.73928\n"
    "            234.00    1.25386     0.76989\n"
    "            268.00    1.24816     0.80050\n"
    "            302.00    1.24283     0.83108 /\n"
    "\n"
    "\n"
    "   94.44    216.50    1.27934     0.67686\n"
    "            250.50    1.27304     0.70995\n"
    "            284.50    1.26716     0.74427\n"
    "            318.50    1.26164     0.77857 /\n"
    "/\n"
    "\n"
    "PVTW\n"
    "    277.0      1.038      4.67E-5    0.318       0.0 /\n"
    "    277.0      1.038      4.67E-5    0.318       0.0 /\n"
    "\n"
    "DENSITY\n"
    "      859.5  1033.0    0.854  /\n"
    "      860.04 1033.0    0.853  /\n"
    "\n";

template <class Evaluation>
inline void testAll()
{
    // test the black-oil specific methods of BlackOilFluidSystem. The generic methods
    // for fluid systems are already tested by the generic test for all fluidsystems.

    typedef typename Opm::MathToolbox<Evaluation>::Scalar Scalar;
    typedef Opm::FluidSystems::BlackOil<double> FluidSystem;

    static constexpr int numPhases = FluidSystem::numPhases;

    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr int gasCompIdx = FluidSystem::gasCompIdx;
    static constexpr int oilCompIdx = FluidSystem::oilCompIdx;
    static constexpr int waterCompIdx = FluidSystem::waterCompIdx;

    Opm::Parser parser;
    Opm::ParseContext parseContext;

    auto deck = parser.parseString(deckString1, parseContext);
    Opm::EclipseState eclState(deck, parseContext);

    FluidSystem::initFromDeck(deck, eclState);

    // create a parameter cache
    typedef typename FluidSystem::template ParameterCache<Scalar> ParamCache;
    ParamCache paramCache(/*maxOilSat=*/0.5, /*regionIdx=*/1);
    if (paramCache.regionIndex() != 1)
        std::abort();

    if (Opm::abs(paramCache.maxOilSat() - 0.5) > 1e-10)
        std::abort();

    if (Opm::abs(FluidSystem::reservoirTemperature() - (273.15 + 15.555)) > 1e-10)
        std::abort();

    if (!FluidSystem::enableDissolvedGas())
        std::abort();

    if (!FluidSystem::enableVaporizedOil())
        std::abort();

    if (FluidSystem::numRegions() != 2)
        std::abort();

    if (FluidSystem::numActivePhases() != 3)
        std::abort();

    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
        if (!FluidSystem::phaseIsActive(phaseIdx))
            std::abort();

    if (FluidSystem::solventComponentIndex(oilPhaseIdx) != oilCompIdx)
        std::abort();
    if (FluidSystem::solventComponentIndex(gasPhaseIdx) != gasCompIdx)
        std::abort();
    if (FluidSystem::solventComponentIndex(waterPhaseIdx) != waterCompIdx)
        std::abort();

    if (FluidSystem::soluteComponentIndex(oilPhaseIdx) != gasCompIdx)
        std::abort();
    if (FluidSystem::soluteComponentIndex(gasPhaseIdx) != oilCompIdx)
        std::abort();

    if (std::abs(FluidSystem::referenceDensity(oilPhaseIdx, /*regionIdx=*/1) - 860.04) > 1e-10)
        std::abort();
    if (std::abs(FluidSystem::referenceDensity(gasPhaseIdx, /*regionIdx=*/1) - 0.853) > 1e-10)
        std::abort();
    if (std::abs(FluidSystem::referenceDensity(waterPhaseIdx, /*regionIdx=*/1) - 1033) > 1e-10)
        std::abort();

    Opm::BlackOilFluidState<Scalar, FluidSystem> fluidState;
    Opm::Valgrind::SetUndefined(fluidState);

    static const Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    unsigned regionIdx = paramCache.regionIndex();
    for (unsigned i = 0; i < 1000; ++i) {
        Scalar p = Scalar(i)/1000*350e5 + 100e5;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState.setPressure(phaseIdx, p);
            fluidState.setSaturation(phaseIdx, 1e-3);
        }

        Scalar RsSat = FluidSystem::saturatedDissolutionFactor(fluidState, oilPhaseIdx, regionIdx);
        fluidState.setRs(RsSat);

        Scalar RvSat = FluidSystem::saturatedDissolutionFactor(fluidState, gasPhaseIdx, regionIdx);
        fluidState.setRv(RvSat);

        paramCache.updateAll(fluidState);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // ensure that the black-oil specific variants of the methods to compute
            // thermdynamic properties return the same value as the generic ones and that
            // the generic methods return the same value as the ones for the saturated
            // quantities (we specify the fluid state to be on the saturation line)
            if (Opm::abs(FluidSystem::density(fluidState, paramCache, phaseIdx)
                         - FluidSystem::density(fluidState, phaseIdx, regionIdx)) > eps)
                std::abort();

            if (Opm::abs(FluidSystem::density(fluidState, paramCache, phaseIdx)
                         - FluidSystem::saturatedDensity(fluidState, phaseIdx, regionIdx)) > eps)
                std::abort();

            Scalar b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, regionIdx);
            Scalar bSat = FluidSystem::saturatedInverseFormationVolumeFactor(fluidState, phaseIdx, regionIdx);
            if (Opm::abs(b - bSat) > eps)
                std::abort();

            if (Opm::abs(FluidSystem::viscosity(fluidState, paramCache, phaseIdx)
                         - FluidSystem::viscosity(fluidState, phaseIdx, regionIdx)) > 1e-10)
                std::abort();

            Scalar R = FluidSystem::saturatedDissolutionFactor(fluidState, phaseIdx, regionIdx);
            Scalar R2 = FluidSystem::saturatedDissolutionFactor(fluidState, phaseIdx, regionIdx);
            if (Opm::abs(R - R2) > eps)
                // seems like there is a problem with D2
                std::abort();

            if (phaseIdx != waterPhaseIdx && // water is immiscible and thus there is no saturation pressure
                Opm::abs(FluidSystem::saturationPressure(fluidState, phaseIdx, regionIdx) - p) > eps*p)
                std::abort();
        }

        if (Opm::abs(FluidSystem::bubblePointPressure(fluidState, regionIdx) - p) > eps*p)
            std::abort();

        if (Opm::abs(FluidSystem::dewPointPressure(fluidState, regionIdx) - p) > eps*p)
            std::abort();

    }

    // make sure that the {oil,gas,water}Pvt() methods are available
    const auto& gPvt OPM_UNUSED = FluidSystem::gasPvt();
    const auto& oPvt OPM_UNUSED  = FluidSystem::oilPvt();
    const auto& wPvt OPM_UNUSED = FluidSystem::waterPvt();
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Opm::DenseAd::Evaluation<double, 2> TestEval;

    testAll<double>();
    //testAll<float>();
    testAll<TestEval>();

    return 0;
}
