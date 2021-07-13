// Created by Christopher Kotthoff on 13.07.2021
// Copyright 2019 ETH Zurich


#include <iostream>

#include <fstream>

#include <mpi.h>

#include <cassert>

#include <cstring>

#include <vector>

#include <math.h>

#include <unordered_map>

#include <unordered_set>

#include <unistd.h>

#include <string>

#include <cmath>

#include <limits.h>

#include <string>

#include <distr/distrbasic.h>
#include <dump/hdf.h>
#include <func/init.h>
#include <parse/argparse.h>
#include <parse/vars.h>
#include <util/distr.h>
#include <util/filesystem.h>
#include <util/format.h>
#include <util/vof.h>

#define CONDITION_A                                               \
  (*fccl[ln])[getCellFromIndex(index, ix - 1, iy - 1, iz - 1)] != \
      kClNone&& fccl_x ==                                         \
      (*fccl[ln])[getCellFromIndex(index, ix - 1, iy - 1, iz - 1)]
#define CONDITION_B                                           \
  (*fccl[ln])[getCellFromIndex(index, ix, iy - 1, iz - 1)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix, iy - 1, iz - 1)]
#define CONDITION_C                                               \
  (*fccl[ln])[getCellFromIndex(index, ix + 1, iy - 1, iz - 1)] != \
      kClNone&& fccl_x ==                                         \
      (*fccl[ln])[getCellFromIndex(index, ix + 1, iy - 1, iz - 1)]
#define CONDITION_D                                           \
  (*fccl[ln])[getCellFromIndex(index, ix - 1, iy, iz - 1)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix - 1, iy, iz - 1)]
#define CONDITION_E                                                           \
  (*fccl[ln])[getCellFromIndex(index, ix, iy, iz - 1)] != kClNone&& fccl_x == \
      (*fccl[ln])[getCellFromIndex(index, ix, iy, iz - 1)]
#define CONDITION_F                                           \
  (*fccl[ln])[getCellFromIndex(index, ix + 1, iy, iz - 1)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix + 1, iy, iz - 1)]
#define CONDITION_G                                               \
  (*fccl[ln])[getCellFromIndex(index, ix - 1, iy + 1, iz - 1)] != \
      kClNone&& fccl_x ==                                         \
      (*fccl[ln])[getCellFromIndex(index, ix - 1, iy + 1, iz - 1)]
#define CONDITION_H                                           \
  (*fccl[ln])[getCellFromIndex(index, ix, iy + 1, iz - 1)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix, iy + 1, iz - 1)]
#define CONDITION_I                                               \
  (*fccl[ln])[getCellFromIndex(index, ix + 1, iy + 1, iz - 1)] != \
      kClNone&& fccl_x ==                                         \
      (*fccl[ln])[getCellFromIndex(index, ix + 1, iy + 1, iz - 1)]
#define CONDITION_P                                           \
  (*fccl[ln])[getCellFromIndex(index, ix - 1, iy - 1, iz)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix - 1, iy - 1, iz)]
#define CONDITION_Q                                                           \
  (*fccl[ln])[getCellFromIndex(index, ix, iy - 1, iz)] != kClNone&& fccl_x == \
      (*fccl[ln])[getCellFromIndex(index, ix, iy - 1, iz)]
#define CONDITION_R                                           \
  (*fccl[ln])[getCellFromIndex(index, ix + 1, iy - 1, iz)] != \
      kClNone&& fccl_x ==                                     \
      (*fccl[ln])[getCellFromIndex(index, ix + 1, iy - 1, iz)]
#define CONDITION_S                                                           \
  (*fccl[ln])[getCellFromIndex(index, ix - 1, iy, iz)] != kClNone&& fccl_x == \
      (*fccl[ln])[getCellFromIndex(index, ix - 1, iy, iz)]

#define Assign_A                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix - 1, iy - 1, iz - 1)]
#define Assign_B                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix, iy - 1, iz - 1)]
#define Assign_C                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix + 1, iy - 1, iz - 1)]
#define Assign_D                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix - 1, iy, iz - 1)]
#define Assign_E                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix, iy, iz - 1)]
#define Assign_F                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix + 1, iy, iz - 1)]
#define Assign_G                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix - 1, iy + 1, iz - 1)]
#define Assign_H                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix, iy + 1, iz - 1)]
#define Assign_I                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix + 1, iy + 1, iz - 1)]
#define Assign_P                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix - 1, iy - 1, iz)]
#define Assign_Q                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix, iy - 1, iz)]
#define Assign_R                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix + 1, iy - 1, iz)]
#define Assign_S                \
  merger_array[merger_size++] = \
      fccl_new[ln][getCellFromIndex(index, ix - 1, iy, iz)]

#define A fccl_new[ln][getCellFromIndex(index, ix - 1, iy - 1, iz - 1)]
#define B fccl_new[ln][getCellFromIndex(index, ix, iy - 1, iz - 1)]
#define C fccl_new[ln][getCellFromIndex(index, ix + 1, iy - 1, iz - 1)]
#define D fccl_new[ln][getCellFromIndex(index, ix - 1, iy, iz - 1)]
#define E fccl_new[ln][getCellFromIndex(index, ix, iy, iz - 1)]
#define F fccl_new[ln][getCellFromIndex(index, ix + 1, iy, iz - 1)]
#define G fccl_new[ln][getCellFromIndex(index, ix - 1, iy + 1, iz - 1)]
#define H fccl_new[ln][getCellFromIndex(index, ix, iy + 1, iz - 1)]
#define I fccl_new[ln][getCellFromIndex(index, ix + 1, iy + 1, iz - 1)]
#define P fccl_new[ln][getCellFromIndex(index, ix - 1, iy - 1, iz)]
#define Q fccl_new[ln][getCellFromIndex(index, ix, iy - 1, iz)]
#define R fccl_new[ln][getCellFromIndex(index, ix + 1, iy - 1, iz)]
#define S fccl_new[ln][getCellFromIndex(index, ix - 1, iy, iz)]

#define MERGE2(first, second)          \
  merger_array[merger_size++] = first; \
  merger_array[merger_size++] = second;

#define MERGE3(first, second, third)    \
  merger_array[merger_size++] = first;  \
  merger_array[merger_size++] = second; \
  merger_array[merger_size++] = third;

#define MERGE4(first, second, third, fourth) \
  merger_array[merger_size++] = first;       \
  merger_array[merger_size++] = second;      \
  merger_array[merger_size++] = third;       \
  merger_array[merger_size++] = fourth;

#define CREATE_3D_DECISION_TREE(                                               \
    useA, useB, useC, useD, useE, useF, useG, useH, useI, useP, useQ, useR,    \
    useS)                                                                      \
  for (auto l : layers) {                                                      \
    auto& fccl_x = (*fccl[l])[getCellFromIndex(index, ix, iy, iz)];            \
    if (DEBUG_CONDITION2 &&                                                    \
        (ix - start[0] > 15 || iy - start[1] > 15 || iz - start[2] > 15))      \
      std::cout << "uh uh out of range... " << std::endl;                      \
    if (DEBUG_CONDITION2 && ix - start[0] == 15 && iy - start[1] == 5 &&       \
        iz - start[2] == 10)                                                   \
      std::cout << "should come now, value of fccl_x is: " << fccl_x           \
                << std::endl;                                                  \
    if (fccl_x != kClNone) {                                                   \
      merger_size = 0;                                                         \
      for (auto ln : layers) {                                                 \
        if (useE && CONDITION_E) {                                             \
          if (DEBUG_CONDITION2 && ix - start[0] == 15 && iy - start[1] == 5 && \
              iz - start[2] == 10)                                             \
            std::cout << "xyz1" << std::endl;                                  \
          Assign_E;                                                            \
        } else {                                                               \
          if (useQ && CONDITION_Q) {                                           \
            if (useH && CONDITION_H) {                                         \
              if (useD && CONDITION_D) {                                       \
                if (DEBUG_CONDITION2 && ix - start[0] == 15 &&                 \
                    iy - start[1] == 5 && iz - start[2] == 10)                 \
                  std::cout << "xyz2" << std::endl;                            \
                Assign_Q;                                                      \
              } else {                                                         \
                if (useF && CONDITION_F) {                                     \
                  if (DEBUG_CONDITION2 && ix - start[0] == 15 &&               \
                      iy - start[1] == 5 && iz - start[2] == 10)               \
                    std::cout << "xyz3" << std::endl;                          \
                  Assign_Q;                                                    \
                } else {                                                       \
                  if (DEBUG_CONDITION2 && ix - start[0] == 15 &&               \
                      iy - start[1] == 5 && iz - start[2] == 10)               \
                    std::cout << "xyz4" << std::endl;                          \
                  MERGE2(H, Q)                                                 \
                }                                                              \
              }                                                                \
            } else {                                                           \
              if (useG && CONDITION_G) {                                       \
                if (useD && CONDITION_D) {                                     \
                  if (useI && CONDITION_I) {                                   \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz5" << std::endl;                        \
                    MERGE2(Q, I)                                               \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz6" << std::endl;                        \
                    Assign_Q;                                                  \
                  }                                                            \
                } else {                                                       \
                  if (useI && CONDITION_I) {                                   \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz7" << std::endl;                        \
                    MERGE3(G, I, Q)                                            \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz8" << std::endl;                        \
                    MERGE2(G, Q)                                               \
                  }                                                            \
                }                                                              \
              } else {                                                         \
                if (useI && CONDITION_I) {                                     \
                  if (useF && CONDITION_F) {                                   \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz9" << std::endl;                        \
                    Assign_Q;                                                  \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz10" << std::endl;                       \
                    MERGE2(I, Q)                                               \
                  }                                                            \
                } else {                                                       \
                  if (DEBUG_CONDITION2 && ix - start[0] == 15 &&               \
                      iy - start[1] == 5 && iz - start[2] == 10)               \
                    std::cout << "xyz11" << std::endl;                         \
                  Assign_Q;                                                    \
                }                                                              \
              }                                                                \
            }                                                                  \
          } else {                                                             \
            if (useB && CONDITION_B) {                                         \
              if (useH && CONDITION_H) {                                       \
                if (useD && CONDITION_D) {                                     \
                  if (DEBUG_CONDITION2 && ix - start[0] == 15 &&               \
                      iy - start[1] == 5 && iz - start[2] == 10)               \
                    std::cout << "xyz12" << std::endl;                         \
                  Assign_B;                                                    \
                                                                               \
                } else {                                                       \
                  if (useF && CONDITION_F) {                                   \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz13" << std::endl;                       \
                    Assign_B;                                                  \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz14" << std::endl;                       \
                    MERGE2(B, H)                                               \
                  }                                                            \
                }                                                              \
              } else {                                                         \
                if (useG && CONDITION_G) {                                     \
                  if (useD && CONDITION_D) {                                   \
                    if (useI && CONDITION_I) {                                 \
                      if (useF && CONDITION_F) {                               \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz15" << std::endl;                   \
                        Assign_B;                                              \
                      } else {                                                 \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz16" << std::endl;                   \
                        MERGE2(B, I)                                           \
                      }                                                        \
                    } else {                                                   \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz17" << std::endl;                     \
                      Assign_B;                                                \
                    }                                                          \
                  } else {                                                     \
                    if (useI && CONDITION_I) {                                 \
                      if (useF && CONDITION_F) {                               \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz18" << std::endl;                   \
                        MERGE2(B, G)                                           \
                      } else {                                                 \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz19" << std::endl;                   \
                        MERGE3(B, G, I)                                        \
                      }                                                        \
                    } else {                                                   \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz20" << std::endl;                     \
                      MERGE2(B, G)                                             \
                    }                                                          \
                  }                                                            \
                } else {                                                       \
                  if (useI && CONDITION_I) {                                   \
                    if (useF && CONDITION_F) {                                 \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz21" << std::endl;                     \
                      Assign_B;                                                \
                    } else {                                                   \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz22" << std::endl;                     \
                      MERGE2(B, I)                                             \
                    }                                                          \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz23" << std::endl;                       \
                    Assign_B;                                                  \
                  }                                                            \
                }                                                              \
              }                                                                \
            } else {                                                           \
              if (useD && CONDITION_D) {                                       \
                if (useF && CONDITION_F) {                                     \
                  if (useH && CONDITION_H) {                                   \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz24" << std::endl;                       \
                    Assign_D;                                                  \
                                                                               \
                  } else {                                                     \
                    if (DEBUG_CONDITION2 && ix - start[0] == 15 &&             \
                        iy - start[1] == 5 && iz - start[2] == 10)             \
                      std::cout << "xyz25" << std::endl;                       \
                    MERGE2(D, F)                                               \
                  }                                                            \
                } else {                                                       \
                  if (useI && CONDITION_I) {                                   \
                    if (useH && CONDITION_H) {                                 \
                      if (useC && CONDITION_C) {                               \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz26" << std::endl;                   \
                        MERGE2(C, D)                                           \
                      } else {                                                 \
                        if (useR && CONDITION_R) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz27" << std::endl;                 \
                          MERGE2(D, R)                                         \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz28" << std::endl;                 \
                          Assign_D;                                            \
                        }                                                      \
                      }                                                        \
                    } else {                                                   \
                      if (useC && CONDITION_C) {                               \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz29" << std::endl;                   \
                        MERGE3(C, D, I)                                        \
                      } else {                                                 \
                        if (useR && CONDITION_R) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz30" << std::endl;                 \
                          MERGE3(D, I, R)                                      \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz31" << std::endl;                 \
                          MERGE2(D, I)                                         \
                        }                                                      \
                      }                                                        \
                    }                                                          \
                  } else {                                                     \
                    if (useC && CONDITION_C) {                                 \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz32" << std::endl;                     \
                      MERGE2(C, D)                                             \
                    } else {                                                   \
                      if (useR && CONDITION_R) {                               \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz33" << std::endl;                   \
                        MERGE2(D, R)                                           \
                      } else {                                                 \
                        if (DEBUG_CONDITION2 && ix - start[0] == 15 &&         \
                            iy - start[1] == 5 && iz - start[2] == 10)         \
                          std::cout << "xyz34" << std::endl;                   \
                        Assign_D;                                              \
                      }                                                        \
                    }                                                          \
                  }                                                            \
                }                                                              \
              } else {                                                         \
                if (useH && CONDITION_H) {                                     \
                  if (useA && CONDITION_A) {                                   \
                    if (useS && CONDITION_S) {                                 \
                      if (useC && CONDITION_C) {                               \
                        if (useF && CONDITION_F) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz35" << std::endl;                 \
                          Assign_H;                                            \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz36" << std::endl;                 \
                          MERGE2(C, H)                                         \
                        }                                                      \
                      } else {                                                 \
                        if (useR && CONDITION_R) {                             \
                          if (useF && CONDITION_F) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz37" << std::endl;               \
                            Assign_H;                                          \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz38" << std::endl;               \
                            MERGE2(H, R)                                       \
                          }                                                    \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz39" << std::endl;                 \
                          Assign_H;                                            \
                        }                                                      \
                      }                                                        \
                    } else {                                                   \
                      if (useC && CONDITION_C) {                               \
                        if (useF && CONDITION_F) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz40" << std::endl;                 \
                          MERGE2(A, H)                                         \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz41" << std::endl;                 \
                          MERGE3(A, C, H)                                      \
                        }                                                      \
                      } else {                                                 \
                        if (useR && CONDITION_R) {                             \
                          if (useF && CONDITION_F) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz42" << std::endl;               \
                            MERGE2(A, H)                                       \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz43" << std::endl;               \
                            MERGE3(A, H, R)                                    \
                          }                                                    \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz44" << std::endl;                 \
                          MERGE2(A, H)                                         \
                        }                                                      \
                      }                                                        \
                    }                                                          \
                  } else {                                                     \
                    if (useP && CONDITION_P) {                                 \
                      if (useS && CONDITION_S) {                               \
                        if (useC && CONDITION_C) {                             \
                          if (useF && CONDITION_F) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz45" << std::endl;               \
                            Assign_H;                                          \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz46" << std::endl;               \
                            MERGE2(C, H)                                       \
                          }                                                    \
                        } else {                                               \
                          if (useR && CONDITION_R) {                           \
                            if (useF && CONDITION_F) {                         \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz47" << std::endl;             \
                              Assign_H;                                        \
                            } else {                                           \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz48" << std::endl;             \
                              MERGE2(H, R)                                     \
                            }                                                  \
                                                                               \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz49" << std::endl;               \
                            Assign_H;                                          \
                          }                                                    \
                        }                                                      \
                      } else {                                                 \
                        if (useC && CONDITION_C) {                             \
                          if (useF && CONDITION_F) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz50" << std::endl;               \
                            MERGE2(H, P)                                       \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz51" << std::endl;               \
                            MERGE3(C, H, P)                                    \
                          }                                                    \
                        } else {                                               \
                          if (useR && CONDITION_R) {                           \
                            if (useF && CONDITION_F) {                         \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz52" << std::endl;             \
                              MERGE2(H, P)                                     \
                            } else {                                           \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz53" << std::endl;             \
                              MERGE3(H, P, R)                                  \
                            }                                                  \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz54" << std::endl;               \
                            MERGE2(H, P)                                       \
                          }                                                    \
                        }                                                      \
                      }                                                        \
                    } else {                                                   \
                      if (useC && CONDITION_C) {                               \
                        if (useF && CONDITION_F) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz55" << std::endl;                 \
                          Assign_H;                                            \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz56" << std::endl;                 \
                          MERGE2(C, H)                                         \
                        }                                                      \
                      } else {                                                 \
                        if (useR && CONDITION_R) {                             \
                          if (useF && CONDITION_F) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz57" << std::endl;               \
                            Assign_H;                                          \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz58" << std::endl;               \
                            MERGE2(H, R)                                       \
                          }                                                    \
                        } else {                                               \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz59" << std::endl;                 \
                          Assign_H;                                            \
                        }                                                      \
                      }                                                        \
                    }                                                          \
                  }                                                            \
                } else {                                                       \
                  if (useF && CONDITION_F) {                                   \
                    if (useS && CONDITION_S) {                                 \
                      if (DEBUG_CONDITION2 && ix - start[0] == 15 &&           \
                          iy - start[1] == 5 && iz - start[2] == 10)           \
                        std::cout << "xyz60" << std::endl;                     \
                      MERGE2(F, S)                                             \
                    } else {                                                   \
                      if (useG && CONDITION_G) {                               \
                        if (useA && CONDITION_A) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz61" << std::endl;                 \
                          MERGE3(A, F, G)                                      \
                        } else {                                               \
                          if (useP && CONDITION_P) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz62" << std::endl;               \
                            MERGE3(F, G, P)                                    \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz63" << std::endl;               \
                            MERGE2(F, G)                                       \
                          }                                                    \
                        }                                                      \
                      } else {                                                 \
                        if (useA && CONDITION_A) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz64" << std::endl;                 \
                          MERGE2(A, F)                                         \
                        } else {                                               \
                          if (useP && CONDITION_P) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz65" << std::endl;               \
                            MERGE2(F, P)                                       \
                                                                               \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz66" << std::endl;               \
                            Assign_F;                                          \
                          }                                                    \
                        }                                                      \
                      }                                                        \
                    }                                                          \
                  } else {                                                     \
                    if (useS && CONDITION_S) {                                 \
                      if (useI && CONDITION_I) {                               \
                        if (useC && CONDITION_C) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz67" << std::endl;                 \
                          MERGE3(C, I, S)                                      \
                        } else {                                               \
                          if (useR && CONDITION_R) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz68" << std::endl;               \
                            MERGE3(I, R, S)                                    \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz69" << std::endl;               \
                            MERGE2(I, S)                                       \
                          }                                                    \
                        }                                                      \
                      } else {                                                 \
                        if (useC && CONDITION_C) {                             \
                          if (DEBUG_CONDITION2 && ix - start[0] == 15 &&       \
                              iy - start[1] == 5 && iz - start[2] == 10)       \
                            std::cout << "xyz70" << std::endl;                 \
                          MERGE2(C, S)                                         \
                        } else {                                               \
                          if (useR && CONDITION_R) {                           \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz71" << std::endl;               \
                            MERGE2(R, S)                                       \
                          } else {                                             \
                            if (DEBUG_CONDITION2 && ix - start[0] == 15 &&     \
                                iy - start[1] == 5 && iz - start[2] == 10)     \
                              std::cout << "xyz72" << std::endl;               \
                            Assign_S;                                          \
                          }                                                    \
                        }                                                      \
                      }                                                        \
                    } else {                                                   \
                      if (useA && CONDITION_A) {                               \
                        if (useC && CONDITION_C) {                             \
                          if (useG && CONDITION_G) {                           \
                            if (useI && CONDITION_I) {                         \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz73" << std::endl;             \
                              MERGE4(A, C, G, I)                               \
                            } else {                                           \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz74" << std::endl;             \
                              MERGE3(A, C, G)                                  \
                            }                                                  \
                          } else {                                             \
                            if (useI && CONDITION_I) {                         \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz75" << std::endl;             \
                              MERGE3(A, C, I)                                  \
                            } else {                                           \
                              if (DEBUG_CONDITION2 && ix - start[0] == 15 &&   \
                                  iy - start[1] == 5 && iz - start[2] == 10)   \
                                std::cout << "xyz76" << std::endl;             \
                              MERGE2(A, C)                                     \
                            }                                                  \
                          }                                                    \
                        } else {                                               \
                          if (useR && CONDITION_R) {                           \
                            if (useG && CONDITION_G) {                         \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz77" << std::endl;           \
                                MERGE4(A, G, I, R)                             \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz78" << std::endl;           \
                                MERGE3(A, G, R)                                \
                              }                                                \
                            } else {                                           \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz79" << std::endl;           \
                                MERGE3(A, I, R)                                \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz80" << std::endl;           \
                                MERGE2(A, R)                                   \
                              }                                                \
                            }                                                  \
                          } else {                                             \
                            if (useG && CONDITION_G) {                         \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz81" << std::endl;           \
                                MERGE3(A, G, I)                                \
                                                                               \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz82" << std::endl;           \
                                MERGE2(A, G)                                   \
                              }                                                \
                            } else {                                           \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz83" << std::endl;           \
                                MERGE2(A, I)                                   \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz84" << std::endl;           \
                                Assign_A;                                      \
                              }                                                \
                            }                                                  \
                          }                                                    \
                        }                                                      \
                      } else {                                                 \
                        if (useC && CONDITION_C) {                             \
                          if (useP && CONDITION_P) {                           \
                            if (useG && CONDITION_G) {                         \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz85" << std::endl;           \
                                MERGE4(C, G, I, P)                             \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz86" << std::endl;           \
                                MERGE3(C, G, P)                                \
                              }                                                \
                            } else {                                           \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz87" << std::endl;           \
                                MERGE3(C, I, P)                                \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz88" << std::endl;           \
                                MERGE2(C, P)                                   \
                              }                                                \
                            }                                                  \
                          } else {                                             \
                            if (useG && CONDITION_G) {                         \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz89" << std::endl;           \
                                MERGE3(C, G, I)                                \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz90" << std::endl;           \
                                MERGE2(C, G)                                   \
                              }                                                \
                            } else {                                           \
                              if (useI && CONDITION_I) {                       \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz91" << std::endl;           \
                                MERGE2(C, I)                                   \
                              } else {                                         \
                                if (DEBUG_CONDITION2 && ix - start[0] == 15 && \
                                    iy - start[1] == 5 && iz - start[2] == 10) \
                                  std::cout << "xyz92" << std::endl;           \
                                Assign_C;                                      \
                              }                                                \
                            }                                                  \
                          }                                                    \
                        } else {                                               \
                          if (useG && CONDITION_G) {                           \
                            if (useI && CONDITION_I) {                         \
                              if (useP && CONDITION_P) {                       \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz93" << std::endl;         \
                                  MERGE4(G, I, P, R)                           \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz94" << std::endl;         \
                                  MERGE3(G, I, P)                              \
                                }                                              \
                              } else {                                         \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz95" << std::endl;         \
                                  MERGE3(G, I, R)                              \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz96" << std::endl;         \
                                  MERGE2(G, I)                                 \
                                }                                              \
                              }                                                \
                            } else {                                           \
                              if (useP && CONDITION_P) {                       \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz97" << std::endl;         \
                                  MERGE3(G, P, R)                              \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz98" << std::endl;         \
                                  MERGE2(G, P)                                 \
                                }                                              \
                              } else {                                         \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz99" << std::endl;         \
                                  MERGE2(G, R)                                 \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz100" << std::endl;        \
                                  Assign_G;                                    \
                                }                                              \
                              }                                                \
                            }                                                  \
                          } else {                                             \
                            if (useI && CONDITION_I) {                         \
                              if (useP && CONDITION_P) {                       \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz101" << std::endl;        \
                                  MERGE3(I, P, R)                              \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz102" << std::endl;        \
                                  MERGE2(I, P)                                 \
                                }                                              \
                              } else {                                         \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz103" << std::endl;        \
                                  MERGE2(I, R)                                 \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz104" << std::endl;        \
                                  Assign_I;                                    \
                                }                                              \
                              }                                                \
                            } else {                                           \
                              if (useP && CONDITION_P) {                       \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz105" << std::endl;        \
                                  MERGE2(P, R)                                 \
                                } else {                                       \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz106" << std::endl;        \
                                  Assign_P;                                    \
                                }                                              \
                              } else {                                         \
                                if (useR && CONDITION_R) {                     \
                                  if (DEBUG_CONDITION2 &&                      \
                                      ix - start[0] == 15 &&                   \
                                      iy - start[1] == 5 &&                    \
                                      iz - start[2] == 10)                     \
                                    std::cout << "xyz107" << std::endl;        \
                                  Assign_R;                                    \
                                } else {                                       \
                                  /*NEW_LABEL_ACTION*/                         \
                                }                                              \
                              }                                                \
                            }                                                  \
                          }                                                    \
                        }                                                      \
                      }                                                        \
                    }                                                          \
                  }                                                            \
                }                                                              \
              }                                                                \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      if (merger_size > 1) {                                                   \
        if (DEBUG_CONDITION2 && merger_array[0] == kClNone)                    \
          std::cout << "Attention1, wrongly assigning kClNone to a cell"       \
                    << std::endl;                                              \
        fccl_new[l][getCellFromIndex(index, ix, iy, iz)] = merger_array[0];    \
        MergeLabels(merger_array, merger_size, union_find_array);              \
      } else if (merger_size == 1) {                                           \
        if (DEBUG_CONDITION2 && merger_array[0] == kClNone)                    \
          std::cout << "Attention2, wrongly assigning kClNone to a cell"       \
                    << std::endl;                                              \
        fccl_new[l][getCellFromIndex(index, ix, iy, iz)] = merger_array[0];    \
      } else {                                                                 \
        fccl_new[l][getCellFromIndex(index, ix, iy, iz)] = ++color;            \
        union_find_array.push_back(color);                                     \
      }                                                                        \
    }                                                                          \
    if (fccl_x == kClNone &&                                                   \
        fccl_new[l][getCellFromIndex(index, ix, iy, iz)] != fccl_x &&          \
        DEBUG_CONDITION2)                                                      \
      std::cout << ix - start[0] << " " << iy - start[1] << " "                \
                << iz - start[2]                                               \
                << " mistake1, fccl_new is not kClNone but fccl is"            \
                << std::endl;                                                  \
    else if (                                                                  \
        fccl_x != kClNone &&                                                   \
        fccl_new[l][getCellFromIndex(index, ix, iy, iz)] == kClNone &&         \
        DEBUG_CONDITION2) {                                                    \
      std::cout << ix - start[0] << " " << iy - start[1] << " "                \
                << iz - start[2]                                               \
                << " mistake2, fccl_new is kClNone but fccl has value. fccl:"  \
                << std::endl;                                                  \
                                                                               \
      for (int ztemp = (iz == start[2] ? iz : iz - 1); ztemp <= iz; ztemp++) { \
        for (int ytemp = start[1]; ytemp <= (ztemp == iz ? iy : end[1] - 1);   \
             ytemp++) {                                                        \
          for (int xtemp = start[0];                                           \
               xtemp <= (ztemp == iz && ytemp == iy ? ix : end[0] - 1);        \
               xtemp++) {                                                      \
            std::cout                                                          \
                << (*fccl[l])[getCellFromIndex(index, xtemp, ytemp, ztemp)]    \
                << " ";                                                        \
          }                                                                    \
          std::cout << std::endl;                                              \
        }                                                                      \
        std::cout << "----------------------------------------" << std::endl;  \
      }                                                                        \
      std::cout << "---------------fccl_new:----------------" << std::endl;    \
      for (int ztemp = (iz == start[2] ? iz : iz - 1); ztemp <= iz; ztemp++) { \
        for (int ytemp = start[1]; ytemp <= (ztemp == iz ? iy : end[1] - 1);   \
             ytemp++) {                                                        \
          for (int xtemp = start[0];                                           \
               xtemp <= (ztemp == iz && ytemp == iy ? ix : end[0] - 1);        \
               xtemp++) {                                                      \
            std::cout                                                          \
                << fccl_new[l][getCellFromIndex(index, xtemp, ytemp, ztemp)]   \
                << " ";                                                        \
          }                                                                    \
          std::cout << std::endl;                                              \
        }                                                                      \
        std::cout << "----------------------------------------" << std::endl;  \
      }                                                                        \
    }                                                                          \
  }

#define CAG_CONSTRUCTION(                                                      \
    equivOfxy_1z_1, equivOfxy_1z, equivOfxy_1z1, equivOfxyz_1, equivOfxyz1,    \
    equivOfxy1z_1, equivOfxy1z, equivOfxy1z1, direction, slowDimVariable,      \
    slowDimIndex, fastDimVariable, fastDimIndex)                               \
  for (auto l : layers) {                                                      \
    { /*case slowDimVariable == 0*/                                            \
      int slowDimVariable = start[slowDimIndex];                               \
      { /*case slowDimVariable == 0, fastDimVariable == 0*/                    \
        int fastDimVariable = start[fastDimIndex];                             \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z1]]         \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      for (int fastDimVariable = start[fastDimIndex] + 1;                      \
           fastDimVariable < end[fastDimIndex] - 1;                            \
           ++fastDimVariable) { /*case slowDimVariable == 0, fastDimVariable== \
                                   middle*/                                    \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z1]]         \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
                                                                               \
      { /*case slowDimVariable == 0, fastDimVariable == last*/                 \
        int fastDimVariable = end[fastDimIndex] - 1;                           \
                                                                               \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
                                                                               \
    for (int slowDimVariable = start[slowDimIndex] + 1;                        \
         slowDimVariable < end[slowDimIndex] - 1;                              \
         ++slowDimVariable) { /*case slowDimVariable == middle*/               \
      { /*case slowDimVariable == middle, fastDimVariable == 0*/               \
        int fastDimVariable = start[fastDimIndex];                             \
                                                                               \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z_1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z_1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z_1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z1]]         \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      for (int fastDimVariable = start[fastDimIndex] + 1;                      \
           fastDimVariable < end[fastDimIndex] - 1;                            \
           ++fastDimVariable) { /*case slowDimVariable ==                      \
                                   middle,fastDimVariable == middle*/          \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z_1] != kClNone &&  \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z_1])            \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z_1]]       \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z_1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z_1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z_1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z1]]         \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
                                                                               \
      { /*case slowDimVariable == middle, fastDimVariable == last*/            \
        int fastDimVariable = end[fastDimIndex] - 1;                           \
                                                                               \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z_1] != kClNone &&  \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z_1])            \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z_1]]       \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz1] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz1])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz1]]          \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
                                                                               \
    { /*case slowDimVariable == last*/                                         \
      int slowDimVariable = end[slowDimIndex] - 1;                             \
      { /*case slowDimVariable == last, fastDimVariable == 0*/                 \
        int fastDimVariable = start[fastDimIndex];                             \
                                                                               \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z_1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z_1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z_1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      for (int fastDimVariable = start[fastDimIndex] + 1;                      \
           fastDimVariable < end[fastDimIndex] - 1;                            \
           ++fastDimVariable) { /*case slowDimVariable == last,                \
                                   fastDimVariable == middle*/                 \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z_1] != kClNone &&  \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z_1])            \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z_1]]       \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z_1] != kClNone &&   \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z_1])             \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z_1]]        \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy1z] != kClNone &&     \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy1z])               \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy1z]]          \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      { /*case slowDimVariable == last, fastDimVariable == last*/              \
        int fastDimVariable = end[fastDimIndex] - 1;                           \
                                                                               \
        auto& border_cell_value = fccl_new[l][getCellFromIndex direction];     \
        \ 
 auto& old_border_cell_value = (*fccl[l])[getCellFromIndex direction];         \
        if (border_cell_value != kClNone) {                                    \
          for (auto ln : layers) {                                             \
            if (fccl_new[ln][getCellFromIndex(index, ix, iy, iz)] !=           \
                    kClNone &&                                                 \
                old_border_cell_value ==                                       \
                    (*fccl[ln])[getCellFromIndex(index, ix, iy, iz)]) {        \
              local_cag[fccl_new[ln][getCellFromIndex(index, ix, iy, iz)]]     \
                  .tryAddEdge(border_cell_value);                              \
            } else {                                                           \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z_1] != kClNone &&  \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z_1])            \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z_1]]       \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxyz_1] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxyz_1])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxyz_1]]         \
                    .tryAddEdge(border_cell_value);                            \
              if (fccl_new[ln][getCellFromIndex equivOfxy_1z] != kClNone &&    \
                  old_border_cell_value ==                                     \
                      (*fccl[ln])[getCellFromIndex equivOfxy_1z])              \
                local_cag[fccl_new[ln][getCellFromIndex equivOfxy_1z]]         \
                    .tryAddEdge(border_cell_value);                            \
            }                                                                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

#define COORD_TO_BLOCKID(x, y, z)      \
  (x) + (y)*m.flags.global_blocks[0] + \
      (z)*m.flags.global_blocks[0] * m.flags.global_blocks[1]

#define DEBUG_CONDITION blockID == 0
#define DEBUG_CONDITION2 blockID == 0

typedef int CAG_Node;

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using IndexCells = GIndex<IdxCell, 3>;

int kBlockSize = -1;
int kDimensions[] = {-1, -1, -1};

constexpr double kClNone = -1;

int mpi_rank_ = -1;
int mpi_size_ = -1;

enum DomainSide {
  kXpositiv = 0,
  kXnegativ = 1,
  kYpositiv = 2,
  kYnegativ = 3,
  kZpositiv = 4,
  kZnegativ = 5
};

struct CAG_NodeProperties {
  CAG_Node CAG_Node_id;
  CAG_Node pointing_to;
  std::vector<CAG_Node> edges;

  std::string ToString() {
    std::string s = "CAG_Node: " + std::to_string(CAG_Node_id) +
                    " pointing_to: " + std::to_string(pointing_to) + " edges: ";

    for (size_t i = 0; i < edges.size(); i++)
      s += std::to_string(edges[i]) + " ";
    return s;
  }
};

struct CAG_NodesPointerTable {
  CAG_Node CAG_Node_id;
  CAG_Node pointing_to;
};

IdxCell getCellFromIndex(const IndexCells& index, int x, int y, int z) {
  MIdx w(x, y, z);
  return index.GetIdx(w);
}

int Find(int color, int* array) {
  if (array[color - 1] != color)
    array[color - 1] = Find(array[color - 1], array);

  return array[color - 1];
}

void DoUnion(int color1, int color2, int* array) {
  int root1 = Find(color1, array);
  int root2 = Find(color2, array);

  if (root1 < root2)
    array[root2 - 1] = root1;
  else
    array[root1 - 1] = root2;
}

void MakeLookupTable(
    int* array, int* lookupArray, int n, std::vector<CAG_Node>& cclabels,
    int color_offset) {
  int* temp = new int[n];

  for (int i = 0; i < n; i++)
    temp[i] = Find(array[i], array);

  for (int currentRank = 1, ranksGiven = 0;
       ranksGiven < n && currentRank < n + 1; currentRank++) {
    int lowest = INT_MAX;
    for (int j = 0; j < n; j++) {
      int currentVal = temp[j];
      if (currentVal < lowest) {
        lowest = currentVal;
      }
    }
    for (int j = 0; j < n; j++) {
      if (temp[j] == lowest) {
        temp[j] = INT_MAX;
        lookupArray[j] = currentRank + color_offset;
        ranksGiven++;
      }
    }
    cclabels.push_back(currentRank + color_offset);
  }

  delete[] temp;
}

void MergeLabels(
    int merger_array[], int merger_size, std::vector<int>& union_find_array) {
  // if (DEBUG_CONDITION2){
  //   std::cout <<"unionfindsize: "<< union_find_array.size()<<" arraysize: "<<
  //   merger_size << " merging "; for (int i = 0; i < merger_size;i++){
  //     std::cout << merger_array[i] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  for (int i = 0; i < merger_size - 1; i++) {
    if (merger_array[i] != merger_array[i + 1]) {
      DoUnion(merger_array[i], merger_array[i + 1], union_find_array.data());
    }
  }
}

/* void DoTwoPass(
    const GRange<size_t>& layers, const Multi<FieldCell<Scal>*>& fccl,
    Multi<FieldCell<Scal>>& fccl_new, std::vector<int>& union_find_array,
    int start_color, std::vector<CAG_Node>& cclabels, M& m) {
  for (auto l : layers) {
    for (auto c : m.AllCells()) {
      if ((*fccl[l])[c] != kClNone) fccl_new[l][c] = start_color++;
    }
  }

  while (true) {
    bool changed = false;
    for (auto l : layers) {
      for (auto c : m.Cells()) {
        if ((*fccl[l])[c] != kClNone) {
          // Update colors with minimum over neighbours
          for (auto cn : m.Stencil(c)) {
            for (auto ln : layers) {
              if ((*fccl[ln])[cn] == (*fccl[l])[c]) {
                if (fccl_new[ln][cn] < fccl_new[l][c]) {
                  changed = true;
                  fccl_new[l][c] = fccl_new[ln][cn];
                }
              }
            }
          }
        }
      }
    }
    if (!changed) {
      break;
    }
  }

  for (auto l : layers) {
    for (auto c : m.Cells()) {
      if (fccl_new[l][c] != kClNone) {
        if (std::find(cclabels.begin(), cclabels.end(), fccl_new[l][c]) ==
            cclabels.end())
          cclabels.push_back(fccl_new[l][c]);
      }
    }
  }
}
 */
void DoTwoPass(
    const GRange<size_t>& layers, const Multi<FieldCell<Scal>*>& fccl,
    Multi<FieldCell<Scal>>& fccl_new, std::vector<int>& union_find_array,
    int start_color, std::vector<CAG_Node>& cclabels, M& m) {
  int color = 0; // first color given will be 1;

  int blockID = m.GetId();

  int merger_size = 0;
  int merger_array[4 * layers.size()];

  const auto& index = m.GetIndexCells();
  const auto& block = m.GetInBlockCells();
  const MIdx start = block.GetBegin();
  const MIdx size = block.GetSize();
  const MIdx end = block.GetEnd();

  { // z:start
    int iz = start[2];

    { // y:start z:start
      int iy = start[1];

      for (auto l : layers) // x:start y:start z:start
        if ((*fccl[l])[getCellFromIndex(index, start[0], start[1], start[2])] !=
            kClNone) {
          fccl_new[l][getCellFromIndex(index, start[0], start[1], start[2])] =
              ++color;
          union_find_array.push_back(color);
        }

      for (int ix = start[0] + 1; ix < end[0];
           ix++) { // x:middle/end y:start z:start
        CREATE_3D_DECISION_TREE(
            false, false, false, false, false, false, false, false, false,
            false, false, false, true);
      }
    }

    for (int iy = start[1] + 1; iy < end[1]; iy++) { // y:middle/end z:start
      { // x:start y:middle/end z:start
        int ix = start[0];
        CREATE_3D_DECISION_TREE(
            false, false, false, false, false, false, false, false, false,
            false, true, true, false);
      }

      for (int ix = start[0] + 1; ix < end[0] - 1;
           ix++) { // x:middle y:middle/end z:start
        CREATE_3D_DECISION_TREE(
            false, false, false, false, false, false, false, false, false, true,
            true, true, true);
      }

      { // x:end y:middle/end z:start
        int ix = end[0] - 1;
        CREATE_3D_DECISION_TREE(
            false, false, false, false, false, false, false, false, false, true,
            true, false, true);
      }
    }
  }

  for (int iz = start[2] + 1; iz < end[2]; iz++) { // z:middle/end
    { // y:start z:middle/end
      int iy = start[1];

      { // x:start y:start z:middle/end
        int ix = start[0];
        CREATE_3D_DECISION_TREE(
            false, false, false, false, true, true, false, true, true, false,
            false, false, false);
      }

      for (int ix = start[0] + 1; ix < end[0] - 1;
           ix++) { // x:middle y:start z:middle/end
        CREATE_3D_DECISION_TREE(
            false, false, false, true, true, true, true, true, true, false,
            false, false, true);
      }

      { // x:end y:start z:middle/end
        int ix = end[0] - 1;
        CREATE_3D_DECISION_TREE(
            false, false, false, true, true, false, true, true, false, false,
            false, false, true);
      }
    }

    for (int iy = start[1] + 1; iy < end[1] - 1;
         iy++) { // y:middle z:middle/end
      { // x:start y:middle z:middle/end
        int ix = start[0];
        CREATE_3D_DECISION_TREE(
            false, true, true, false, true, true, false, true, true, false,
            true, true, false);
      }

      for (int ix = start[0] + 1; ix < end[0] - 1;
           ix++) { // x:middle y:middle z:middle/end
        CREATE_3D_DECISION_TREE(
            true, true, true, true, true, true, true, true, true, true, true,
            true, true);
      }

      { // x:end y:middle z:middle/end
        int ix = end[0] - 1;
        CREATE_3D_DECISION_TREE(
            true, true, false, true, true, false, true, true, false, true, true,
            false, true);
      }
    }

    { // y:end z:middle/end
      int iy = end[1] - 1;

      { // x:start y:end z:middle/end
        int ix = start[0];
        CREATE_3D_DECISION_TREE(
            false, true, true, false, true, true, false, false, false, false,
            true, true, false);
      }

      for (int ix = start[0] + 1; ix < end[0] - 1;
           ix++) { // x:middle y:end z:middle/end
        CREATE_3D_DECISION_TREE(
            true, true, true, true, true, true, false, false, false, true, true,
            true, true);
      }

      { // x:end y:end z:middle/end
        int ix = end[0] - 1;
        CREATE_3D_DECISION_TREE(
            true, true, false, true, true, false, false, false, false, true,
            true, false, true);
      }
    }
  }

  int* look_up_table = new int[union_find_array.size()];
  // if (mpi_rank_ == 13)
  //		std::cout << "size: " << union_find_array.size() << std::endl;
  MakeLookupTable(
      union_find_array.data(), look_up_table, union_find_array.size(), cclabels,
      start_color);

  // std::cout << std::endl << std::endl << "union_find_array:";
  // for (size_t i = 0; i < union_find_array.size();i++)
  //   std::cout << union_find_array[i] << " ";

  // std::cout << std::endl << std::endl << "look_up_table:";
  // for (size_t i = 0; i < union_find_array.size();i++)
  //   std::cout << look_up_table[i] << " ";
  // if (DEBUG_CONDITION2){
  //   for (auto l : layers) {
  //     for (auto c : m.Cells()) {
  //       if (fccl_new[l][c] != kClNone) {
  //         std::cout << (int) fccl_new[l][c] << " ";
  //       }
  //     }
  //   }
  //   std::cout << std::endl;
  // }

  for (auto l : layers) {
    for (auto c : m.Cells()) {
      if (fccl_new[l][c] != kClNone) {
        fccl_new[l][c] = look_up_table[(int)(fccl_new[l][c]) - 1];
      }
    }
  }

  //  if (DEBUG_CONDITION2){
  //   for (auto l : layers) {
  //     for (auto c : m.Cells()) {
  //       if (fccl_new[l][c] != kClNone) {
  //         std::cout << (int) fccl_new[l][c] << " ";
  //       }
  //     }
  //   }
  //   std::cout << std::endl;
  // }

  delete[] look_up_table;
}

int CagFind(
    int CAG_Node_id, std::unordered_map<CAG_Node, CAG_NodeProperties>& cag) {
  if (cag.at(CAG_Node_id).pointing_to != CAG_Node_id)
    cag[CAG_Node_id].pointing_to =
        CagFind(cag.at(CAG_Node_id).pointing_to, cag);

  return cag.at(CAG_Node_id).pointing_to;
}

void CagDoUnion(
    CAG_Node CAG_Node1, CAG_Node CAG_Node2,
    std::unordered_map<CAG_Node, CAG_NodeProperties>& cag) {
  CAG_Node root1 = CagFind(CAG_Node1, cag);
  CAG_Node root2 = CagFind(CAG_Node2, cag);

  if (root1 == root2) {
    // if (DEBUG_CONDITION) {
    //   std::cout << "  └─nothing to do since roots are the same" << std::endl;
    // }
    return;
  } else if (root1 < root2) {
    cag[root2].pointing_to = root1;

    // if (DEBUG_CONDITION) {
    //   std::cout << "  └─root of CAG_Node " << CAG_Node1 << " is " << root1
    //             << std::endl;
    //   std::cout << "  └─root of CAG_Node " << CAG_Node2 << " is " << root2
    //             << std::endl;
    //   std::cout << "  └─setting pointing_to of root " << root2 << " to "
    //             << root1 << std::endl;
    // }

    // return root1;
  } else {
    cag[root1].pointing_to = root2;

    // if (DEBUG_CONDITION) {
    //   std::cout << "  └─root of CAG_Node " << CAG_Node1 << " is " << root1
    //             << std::endl;
    //   std::cout << "  └─root of CAG_Node " << CAG_Node2 << " is " << root2
    //             << std::endl;
    //   std::cout << "  └─setting pointing_to of root " << root1 << " to "
    //             << root2 << std::endl;
    // }
    // return root2;
  }
}

void RecolorDistributed(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
    Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfc, bool verb,
    bool unionfind, bool reduce, bool grid, M& m) {
  auto sem = m.GetSem("recolor");
  int blockID = m.GetId();
  std::cout << m.GetId() << std::endl;
  if (DEBUG_CONDITION2) {
    std::cout << blockID << ": condition2:-2" << std::endl;
  }
  struct {
    std::vector<CAG_Node> cclabels;
    std::unordered_map<CAG_Node, CAG_NodesPointerTable>
        local_CAG_Nodes_pointer_table;
    std::vector<int> partners;
    std::unordered_map<CAG_Node, CAG_NodeProperties> cag;
    int partners_size;
    int counter;
    int lowest_local_CAG_Node_id;
    int highest_local_CAG_Node_id;
    int local_cc_size;
    Multi<FieldCell<Scal>> fccl_new;
    std::vector<int> compressed_cag;
    std::vector<int> received_compressed_cag;

    int neighbours[6];

    // Pointers to objects from other local blocks:
    std::vector<std::vector<int>*> collected_compressed_cags;
    std::vector<std::vector<int>*> collected_recieved_compressed_cags;

    // only relevant to lead block:
    bool first_reduction;
    std::unordered_map<int, std::vector<int>*> receiver_cag_map;

  } * ctx(sem);

  if (DEBUG_CONDITION2) {
    std::cout << blockID << ": condition2:-1" << std::endl;
  }
  auto& t = *ctx;

  auto& cclabels = ctx->cclabels;
  auto& local_CAG_Nodes_pointer_table = ctx->local_CAG_Nodes_pointer_table;
  auto& partners = ctx->partners;
  auto& cag = ctx->cag;
  auto& partners_size = ctx->partners_size;
  auto& counter = ctx->counter;
  auto& fccl_new = ctx->fccl_new;
  auto& compressed_cag = ctx->compressed_cag;
  auto& received_compressed_cag = ctx->received_compressed_cag;
  auto& collected_compressed_cags = ctx->collected_compressed_cags;
  auto& collected_recieved_compressed_cags =
      ctx->collected_recieved_compressed_cags;
  auto& first_reduction = ctx->first_reduction;
  auto& receiver_cag_map = ctx->receiver_cag_map;
  auto& neighbours = ctx->neighbours;

  auto& lowest_local_CAG_Node_id = ctx->lowest_local_CAG_Node_id;
  auto& highest_local_CAG_Node_id = ctx->highest_local_CAG_Node_id;
  auto& local_cc_size = ctx->local_cc_size;
  if (DEBUG_CONDITION2) {
    std::cout << blockID << ": condition2:0" << std::endl;
  }

  if (m.IsLead()) {
    MPI_Comm_rank(m.GetMpiComm(), &mpi_rank_);
    MPI_Comm_size(m.GetMpiComm(), &mpi_size_);
  }
  if (DEBUG_CONDITION2) {
    std::cout << blockID << ": condition2:1" << std::endl;
  }
  kDimensions[0] = m.GetGlobalSize().data()[0];
  kDimensions[1] = m.GetGlobalSize().data()[1];
  kDimensions[2] = m.GetGlobalSize().data()[2];

  kBlockSize = m.GetInBlockCells().GetSize()[0];

  if (DEBUG_CONDITION2) {
    std::cout << blockID << ": condition2:2" << std::endl;
  }
  int start_color =
      blockID * kBlockSize * kBlockSize * kBlockSize * layers.size();
  if (sem()) {
    // works only with mpi_size == blocks atm
    // assert(
    //     mpi_size_ == m.flags.global_blocks[0] * m.flags.global_blocks[1] *
    //                      m.flags.global_blocks[2]);
    // assert(mpi_rank_ == blockID);
    assert(kDimensions[0] % kBlockSize == 0);
    assert(kDimensions[1] % kBlockSize == 0);
    assert(kDimensions[2] % kBlockSize == 0);
    // if (mpi_size_ != m.flags.global_blocks[0] * m.flags.global_blocks[1] *
    //                      m.flags.global_blocks[2])
    //   std::cout << "ERROR: Amount of ranks should be "
    //             << m.flags.global_blocks[0] * m.flags.global_blocks[1] *
    //                    m.flags.global_blocks[2]
    //             << " but is " << mpi_size_ << "\n";
    // assert(
    //     mpi_size_ == m.flags.global_blocks[0] * m.flags.global_blocks[1] *
    //                      m.flags.global_blocks[2]);
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:3" << std::endl;
    }
    int my_x = blockID % m.flags.global_blocks[0];
    int my_y =
        (blockID % (m.flags.global_blocks[1] * m.flags.global_blocks[0])) /
        m.flags.global_blocks[0];
    int my_z = blockID / (m.flags.global_blocks[1] * m.flags.global_blocks[0]);

    if (DEBUG_CONDITION)
      std::cout << "amount of blocks in xyz: " << m.flags.global_blocks[0]
                << " " << m.flags.global_blocks[1] << " "
                << m.flags.global_blocks[2] << std::endl;
    if (DEBUG_CONDITION)
      std::cout << blockID << ": my block coordinates are " << my_x << " "
                << my_y << " " << my_z << std::endl;
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:4" << std::endl;
    }
    collected_compressed_cags.push_back(
        &compressed_cag); // do once since reference stays the same during all
                          // iterations of the reduction.
    m.GatherToLead(&collected_compressed_cags);
    received_compressed_cag.push_back(blockID);
    collected_recieved_compressed_cags.push_back(
        &received_compressed_cag); // do once since reference stays the same
                                   // during all iterations of the reduction.
    m.GatherToLead(&collected_recieved_compressed_cags);
    first_reduction = true;
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:5" << std::endl;
    }

    {
      if (my_x + 1 < m.flags.global_blocks[0])
        neighbours[kXpositiv] = COORD_TO_BLOCKID(my_x + 1, my_y, my_z);
      else
        neighbours[kXpositiv] =
            COORD_TO_BLOCKID(my_x + 1 - m.flags.global_blocks[0], my_y, my_z);

      if (my_x - 1 >= 0)
        neighbours[kXnegativ] = COORD_TO_BLOCKID(my_x - 1, my_y, my_z);
      else
        neighbours[kXnegativ] =
            COORD_TO_BLOCKID(my_x - 1 + m.flags.global_blocks[0], my_y, my_z);

      if (my_y + 1 < m.flags.global_blocks[1])
        neighbours[kYpositiv] = COORD_TO_BLOCKID(my_x, my_y + 1, my_z);
      else
        neighbours[kYpositiv] =
            COORD_TO_BLOCKID(my_x, my_y + 1 - m.flags.global_blocks[1], my_z);

      if (my_y - 1 >= 0)
        neighbours[kYnegativ] = COORD_TO_BLOCKID(my_x, my_y - 1, my_z);
      else
        neighbours[kYnegativ] =
            COORD_TO_BLOCKID(my_x, my_y - 1 + m.flags.global_blocks[1], my_z);

      if (my_z + 1 < m.flags.global_blocks[2])
        neighbours[kZpositiv] = COORD_TO_BLOCKID(my_x, my_y, my_z + 1);
      else
        neighbours[kZpositiv] =
            COORD_TO_BLOCKID(my_x, my_y, my_z + 1 - m.flags.global_blocks[2]);

      if (my_z - 1 >= 0)
        neighbours[kZnegativ] = COORD_TO_BLOCKID(my_x, my_y, my_z - 1);
      else
        neighbours[kZnegativ] =
            COORD_TO_BLOCKID(my_x, my_y, my_z - 1 + m.flags.global_blocks[2]);
    }
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:6" << std::endl;
    }
  }
  if (sem()) { // two pass and border exchange
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:7" << std::endl;
    }
    std::vector<int> union_find_array;

    fccl_new.Reinit(layers, m, kClNone);

    counter = 0;
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 0 << " " << counter++ << std::endl;

    // if (DEBUG_CONDITION2) {
    //   for (auto c : m.AllCells()) {
    //     bool temp = false;
    //     if (!m.IsInner(c))
    //       for (auto cn : m.Stencil(c))
    //         if (m.IsInner(cn)) temp = true;
    //     if (temp)
    //       for (auto l : layers)
    //         std::cout << (*fccl[l])[c] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:8" << std::endl;
    }

    DoTwoPass(
        layers, fccl, fccl_new, union_find_array, start_color++, cclabels, m);
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 1 << " " << counter++ << std::endl;
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:9" << std::endl;
    }
    // if (DEBUG_CONDITION2) {
    //   for (auto c : m.AllCells()) {
    //     bool temp = false;
    //     if (!m.IsInner(c))
    //       for (auto cn : m.Stencil(c))
    //         if (m.IsInner(cn)) temp = true;
    //     if (temp)
    //       for (auto l : layers)
    //         std::cout << fccl_new[l][c] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:10" << std::endl;
    }
    for (auto l : layers) {
      m.Comm(&fccl_new[l]);
    }
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 2 << " " << counter++ << std::endl;
  }

  if (sem()) { // local cag construction
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:11" << std::endl;
    }
    local_cc_size = cclabels.size();
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 3 << " " << counter++ << std::endl;

    struct Edges {
      CAG_Node CAG_Node_id;
      CAG_Node previous_try_cache;
      std::vector<CAG_Node> edge_to;

      Edges() : previous_try_cache(-1) {}

      void tryAddEdge(CAG_Node color) {
        // if (DEBUG_CONDITION && color == 12386)
        //   std::cout << "HELP, ADDING 12386 rn :(" << std::endl;
        // if already in there, return without adding
        if (color == previous_try_cache) return;
        for (size_t i = 0; i < edge_to.size(); i++)
          if (edge_to[i] == color) {
            previous_try_cache = color;
            return;
          }
        // if not in there, add
        previous_try_cache = color;
        edge_to.push_back(color);
      }
    };
    std::unordered_map<CAG_Node, Edges>
        local_cag; // maybe double instead of cag_node
    for (int i = 0; i < local_cc_size; i++) {
      struct Edges local_node;
      local_node.CAG_Node_id = cclabels[i];

      local_cag[cclabels[i]] = local_node;
    }
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:12" << std::endl;
    }

    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 4 << " " << counter++ << std::endl;

    // if (DEBUG_CONDITION)
    //     for (int i = 0; i < local_cc_size; i++) {
    //       local_cag[cclabels[i]].CAG_Node_id = cclabels[i];
    //       std::cout << local_cag[cclabels[i]].CAG_Node_id << " = " <<
    //       cclabels[i]
    //                  << std::endl;
    //     }
    // enum DomainSide {kXpositiv=0, kXnegativ=1, kYpositiv=2, Ynegativ=3,
    // kZpositiv=4, kZnegativ=5};

    // local cag construction

    const auto& index = m.GetIndexCells();
    const auto& block = m.GetInBlockCells();
    const MIdx start = block.GetBegin();
    const MIdx size = block.GetSize();
    const MIdx end = block.GetEnd();
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 5 << " " << counter++ << std::endl;

    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:13" << std::endl;
    }

    if (DEBUG_CONDITION) {
      std::cout << "am in local_cag, size is: " << local_cag.size()
                << std::endl;
      for (auto it = local_cag.begin(); it != local_cag.end(); ++it) {
        std::cout << "line 1637: cag key " << it->first << " nodeid "
                  << it->second.CAG_Node_id << " edges to: ";

        auto& temp = it->second.edge_to;

        for (int i = 0; i < temp.size(); i++)
          std::cout << temp[i] << " ";
        std::cout << std::endl;
      }
    }

    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:14" << std::endl;
    }

    {
      if (neighbours[kXpositiv] != -1) {
        int ix = end[0] - 1;
        CAG_CONSTRUCTION(
            (index, ix, iy - 1, iz - 1), (index, ix, iy - 1, iz),
            (index, ix, iy - 1, iz + 1), (index, ix, iy, iz - 1),
            (index, ix, iy, iz + 1), (index, ix, iy + 1, iz - 1),
            (index, ix, iy + 1, iz), (index, ix, iy + 1, iz + 1),
            (index, ix + 1, iy, iz), iz, 2, iy, 1);
      }
      if (neighbours[kXnegativ] != -1) {
        int ix = start[0];
        CAG_CONSTRUCTION(
            (index, ix, iy - 1, iz - 1), (index, ix, iy - 1, iz),
            (index, ix, iy - 1, iz + 1), (index, ix, iy, iz - 1),
            (index, ix, iy, iz + 1), (index, ix, iy + 1, iz - 1),
            (index, ix, iy + 1, iz), (index, ix, iy + 1, iz + 1),
            (index, ix - 1, iy, iz), iz, 2, iy, 1);
      }
      if (neighbours[kYpositiv] != -1) {
        int iy = end[1] - 1;
        CAG_CONSTRUCTION(
            (index, ix - 1, iy, iz - 1), (index, ix - 1, iy, iz),
            (index, ix - 1, iy, iz + 1), (index, ix, iy, iz - 1),
            (index, ix, iy, iz + 1), (index, ix + 1, iy, iz - 1),
            (index, ix + 1, iy, iz), (index, ix + 1, iy, iz + 1),
            (index, ix, iy + 1, iz), iz, 2, ix, 0);
      }
      if (neighbours[kYnegativ] != -1) {
        int iy = start[1];
        CAG_CONSTRUCTION(
            (index, ix - 1, iy, iz - 1), (index, ix - 1, iy, iz),
            (index, ix - 1, iy, iz + 1), (index, ix, iy, iz - 1),
            (index, ix, iy, iz + 1), (index, ix + 1, iy, iz - 1),
            (index, ix + 1, iy, iz), (index, ix + 1, iy, iz + 1),
            (index, ix, iy - 1, iz), iz, 2, ix, 0);
      }
      if (neighbours[kZpositiv] != -1) {
        int iz = end[2] - 1;
        CAG_CONSTRUCTION(
            (index, ix - 1, iy - 1, iz), (index, ix - 1, iy, iz),
            (index, ix - 1, iy + 1, iz), (index, ix, iy - 1, iz),
            (index, ix, iy + 1, iz), (index, ix + 1, iy - 1, iz),
            (index, ix + 1, iy, iz), (index, ix + 1, iy + 1, iz),
            (index, ix, iy, iz + 1), iy, 1, ix, 0);
      }
      if (neighbours[kZnegativ] != -1) {
        int iz = start[2];
        CAG_CONSTRUCTION(
            (index, ix - 1, iy - 1, iz), (index, ix - 1, iy, iz),
            (index, ix - 1, iy + 1, iz), (index, ix, iy - 1, iz),
            (index, ix, iy + 1, iz), (index, ix + 1, iy - 1, iz),
            (index, ix + 1, iy, iz), (index, ix + 1, iy + 1, iz),
            (index, ix, iy, iz - 1), iy, 1, ix, 0);
      }
    }

    if (DEBUG_CONDITION) {
      std::cout << "am in local_cag, size is: " << local_cag.size()
                << std::endl;
      for (auto it = local_cag.begin(); it != local_cag.end(); ++it) {
        std::cout << "line 1637: cag key " << it->first << " nodeid "
                  << it->second.CAG_Node_id << " edges to: ";

        auto& temp = it->second.edge_to;

        for (int i = 0; i < temp.size(); i++)
          std::cout << temp[i] << " ";
        std::cout << std::endl;
      }
    }

    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 6 << " " << counter++ << std::endl;

    lowest_local_CAG_Node_id = INT_MAX;
    highest_local_CAG_Node_id = INT_MIN;

    for (int i = 0; i < local_cc_size; i++) { // maybe not correctly
                                              // initialized?
      struct CAG_NodesPointerTable temp;
      temp.CAG_Node_id = cclabels[i];
      temp.pointing_to = cclabels[i];

      local_CAG_Nodes_pointer_table[cclabels[i]] = temp;
      // local_CAG_Nodes_pointer_table[cclabels[i]].CAG_Node_id = cclabels[i];
      // local_CAG_Nodes_pointer_table[cclabels[i]].pointing_to = cclabels[i];
      if (cclabels[i] < lowest_local_CAG_Node_id)
        lowest_local_CAG_Node_id = cclabels[i];
      if (cclabels[i] > highest_local_CAG_Node_id)
        highest_local_CAG_Node_id = cclabels[i];
    }
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:15" << std::endl;
    }
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 7 << " " << counter++ << std::endl;

    // calculating reduction tree

    int blocks = m.flags.global_blocks[0] * m.flags.global_blocks[1] *
                 m.flags.global_blocks[2];

    partners_size = (int)(log(blocks) / log(2.0));
    assert(
        (int)pow(2, partners_size) == blocks &&
        "Amount of blocks must be power of 2 for the reduction tree "
        "(temporarily)");
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 8 << " " << counter++ << std::endl;

    {
      int nxt_distance = 1;
      while (nxt_distance < blocks) {
        int temp = int((blockID) / nxt_distance);
        int skip = 0;

        if (temp % 2 == 1)
          skip = -nxt_distance;
        else
          skip = nxt_distance;

        partners.push_back(blockID + skip);

        nxt_distance *= 2;
      }
    }
    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 9 << " " << counter++ << std::endl;
    if (DEBUG_CONDITION2) {
      std::cout << blockID << ": condition2:16" << std::endl;
    }
    // putting local_cag into cag

    for (auto it = local_cag.begin(); it != local_cag.end(); ++it) {
      if ((it->second).edge_to.size() ==
          0) // dont put it into the cag if there are no outgoing edges
        continue;

      struct CAG_NodeProperties temp;
      temp.CAG_Node_id = (it->second).CAG_Node_id;
      temp.pointing_to = (it->second).CAG_Node_id;

      // std::cout << (it->second).CAG_Node_id << std::endl;
      // std::cout << temp.CAG_Node_id << std::endl;

      for (size_t j = 0; j < (it->second).edge_to.size(); j++) {
        temp.edges.push_back((it->second).edge_to[j]);
        // std::cout << blockID << " putting :" << (it->second).edge_to[j] <<
        // std::endl;
      }

      cag[(it->second).CAG_Node_id] = temp;
    }

    if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 10 << " " << counter++ << std::endl;
  }
  // the actual reduction
  for (int index = 0; index < partners_size; index++) {
    if (sem("reduction")) {
      std::cout << blockID << ": PARTNERS: " << partners_size << std::endl;
      int& partner = partners[index];
      if (DEBUG_CONDITION2) {
        std::cout << blockID << ": condition2:17" << std::endl;
      }
      if (DEBUG_CONDITION) {
        std::cout
            << std::endl
            << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░";
        if (partner >= 10) std::cout << "░";
        std::cout << std::endl
                  << "░░░░░░░░░░░░░░░░ Partner this iteration: " << partner
                  << " ░░░░░░░░░░░░░░░░" << std::endl;
      }
      if (DEBUG_CONDITION2)
        std::cout << blockID << ": " << 11 << " " << counter++ << std::endl;

      if (DEBUG_CONDITION) {
        std::cout << std::endl << "current cag:" << std::endl;
        for (auto it = cag.begin(); it != cag.end(); ++it) {
          std::cout << "└─" << it->second.ToString() << std::endl;
        }
        std::cout << std::endl;
      }
      if (DEBUG_CONDITION2)
        std::cout << blockID << ": " << 12 << " " << counter++ << std::endl;

      // compress cag, to send it over
      compressed_cag.clear();
      compressed_cag.push_back(partner);
      for (auto it = cag.begin(); it != cag.end(); ++it) {
        struct CAG_NodeProperties& temp = it->second;
        compressed_cag.push_back(temp.CAG_Node_id);
        compressed_cag.push_back(temp.pointing_to);

        if (DEBUG_CONDITION2) {
          std::cout << blockID << ": condition2:18" << std::endl;
        }

        if (DEBUG_CONDITION)
          std::cout << "just put into the compressed cag node: "
                    << compressed_cag.back() << std::endl;

        for (size_t j = 0; j < temp.edges.size(); j++)
          if (temp.edges[j] != -1) compressed_cag.push_back(temp.edges[j]);

        compressed_cag.push_back(-1);
      }
      if (DEBUG_CONDITION2)
        std::cout << blockID << ": " << 13 << " " << counter++ << std::endl;

      // send over compressed cag
      // MPI_Request req;
      // MPI_Isend(
      //     compressed_cag.data(), compressed_cag.size(), MPI_INT,
      //     partners[index], 99, m.GetMpiComm(), &req);
      // if (DEBUG_CONDITION2)
      //   std::cout << blockID << ": " << 14 << " " << counter++ << std::endl;
    }
    if (sem("data_exchange") && m.IsLead()) {
      if (DEBUG_CONDITION2) {
        std::cout << blockID << ": condition2:19" << std::endl;
      }
      if (first_reduction) { // is only true for the very first time.
        first_reduction = false;
        for (int i = 0; i < collected_recieved_compressed_cags.size(); i++) {
          receiver_cag_map[(*collected_recieved_compressed_cags[i])[0]] =
              collected_recieved_compressed_cags[i];
        }
      }

      std::vector<std::vector<int>> send_queue(mpi_size_);

      for (auto it = collected_compressed_cags.begin();
           it != collected_compressed_cags.end(); ++it) {
        auto& current_comporessed_cag = **it;
        int rank_to_send_to = m.GetMpiRankFromId(current_comporessed_cag[0]);
        auto& current_send_queue = send_queue[rank_to_send_to];

        current_send_queue.insert(
            current_send_queue.end(), current_comporessed_cag.begin(),
            current_comporessed_cag.end());
        current_send_queue.push_back(-2);
      }
      if (DEBUG_CONDITION2) {
        std::cout << blockID << ": condition2:20" << std::endl;
      }
      std::vector<MPI_Request> requests;

      for (int i = 0; i < send_queue.size(); i++) {
        if (send_queue[i].size() == 0 || i == mpi_rank_) continue;

        requests.emplace_back();
        MPI_Isend(
            send_queue[i].data(), send_queue[i].size(), MPI_INT, i, 969,
            m.GetMpiComm(), &requests.back());
      }
      if (DEBUG_CONDITION2) {
        std::cout << blockID << ": condition2:21" << std::endl;
      }
      for (int i = 0; i < send_queue.size(); i++) {
        if (send_queue[i].size() == 0 || i == mpi_rank_) continue;

        MPI_Status recv_status;
        MPI_Probe(i, 969, m.GetMpiComm(), &recv_status);
        int recv_size = -1;
        MPI_Get_count(&recv_status, MPI_INT, &recv_size);
        std::vector<int> recieved_cag(recv_size);
        MPI_Recv(
            recieved_cag.data(), recv_size, MPI_INT, i, 969, m.GetMpiComm(),
            MPI_STATUS_IGNORE);

        int index = 0;
        while (index < recv_size) {
          int receiver_block = recieved_cag[index++];
          auto& current_received_compressed_cag =
              *(receiver_cag_map.at(receiver_block));

          current_received_compressed_cag.clear();

          while (true) {
            if (recieved_cag[index++] == -2) break;

            current_received_compressed_cag.push_back(recieved_cag[index - 1]);
          }
        }
      }
      if (DEBUG_CONDITION2) {
        std::cout << blockID << ": condition2:22" << std::endl;
      }
      auto& recieved_cag = send_queue[mpi_rank_];
      int index = 0;
      while (index < recieved_cag.size()) {
        int receiver_block = recieved_cag[index++];
        auto& current_received_compressed_cag =
            *(receiver_cag_map.at(receiver_block));

        current_received_compressed_cag.clear();

        while (true) {
          if (recieved_cag[index++] == -2) break;

          current_received_compressed_cag.push_back(recieved_cag[index - 1]);
        }
      }

      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }
    if (sem("receive_and_remainder")) {
      // receiving cag
      // MPI_Status recv_status;
      // MPI_Probe(partners[index], 99, m.GetMpiComm(), &recv_status);
      // int recv_size = -1;
      // MPI_Get_count(&recv_status, MPI_INT, &recv_size);
      // std::vector<int> recieved_cag(recv_size);
      // MPI_Recv(
      //     recieved_cag.data(), recv_size, MPI_INT, partners[index], 99,
      //     m.GetMpiComm(), MPI_STATUS_IGNORE);
      // if (DEBUG_CONDITION2)
      std::cout << blockID << ": " << 15 << " " << counter++ << std::endl;

      // decompressing and unioning received cag
      if (DEBUG_CONDITION) {
        std::cout << "received cag:" << std::endl;
      }
      for (int i = 0; i < received_compressed_cag.size();) {
        struct CAG_NodeProperties temp;
        temp.CAG_Node_id = received_compressed_cag[i++];
        temp.pointing_to = received_compressed_cag[i++];

        while (true) {
          CAG_Node CAG_Node = received_compressed_cag[i++];
          if (CAG_Node == -1) break;
          temp.edges.push_back(CAG_Node);
        }
        if (DEBUG_CONDITION) {
          std::cout << "└─" << temp.ToString() << std::endl;
        }
        cag[temp.CAG_Node_id] = temp;
      }

      if (DEBUG_CONDITION) {
        std::cout << std::endl << "current cag:" << std::endl;
        for (auto it = cag.begin(); it != cag.end(); ++it) {
          std::cout << "└─" << it->second.ToString() << std::endl;
        }
        std::cout << std::endl;
      }
      if (DEBUG_CONDITION2)
        std::cout << blockID << ": " << 16 << " " << counter++ << std::endl;

      // contracting all possible edges
      if (DEBUG_CONDITION) {
        std::cout << std::endl
                  << "contracting all possible edges:" << std::endl;
      }
      std::unordered_set<CAG_Node> has_outgoing_edges_left;
      for (auto it = cag.begin(); it != cag.end(); ++it) {
        struct CAG_NodeProperties& temp = it->second;
        for (size_t i = 0; i < temp.edges.size(); i++) {
          CAG_Node& from = temp.CAG_Node_id;
          CAG_Node& to = temp.edges[i];

          if (to == -1) continue;

          // check if it's an outgoing edge, if yes, abort
          auto it_to = cag.find(to);
          if (it_to == cag.end()) {
            if (DEBUG_CONDITION) {
              std::cout << "└─detected outgoing edge from " << from << " to "
                        << to << std::endl;
            }
            has_outgoing_edges_left.insert(temp.CAG_Node_id);
            continue;
          }

          if (DEBUG_CONDITION) {
            std::cout << "└─unioning " << from << " and " << to << std::endl;
          }
          CagDoUnion(from, to, cag);
          to = -1;
        }
      }

      if (DEBUG_CONDITION) {
        std::cout << std::endl << "current cag:" << std::endl;
        for (auto it = cag.begin(); it != cag.end(); ++it) {
          std::cout << "└─" << it->second.ToString() << std::endl;
        }
        std::cout << std::endl;
      }
      if (DEBUG_CONDITION2)
        std::cout << blockID << ": " << 17 << " " << counter++ << std::endl;

      // remove all components with no outgoing edges
      if (DEBUG_CONDITION) {
        std::cout << std::endl
                  << "remove all components with no outgoing edges:"
                  << std::endl;
      }
      std::unordered_set<CAG_Node> roots_outgoing_edges_left;
      for (auto it = has_outgoing_edges_left.cbegin();
           it != has_outgoing_edges_left.cend(); ++it) {
        roots_outgoing_edges_left.insert(CagFind(*it, cag));
      }
      for (auto it = cag.begin(); it != cag.end(); ++it) {
        CAG_Node& CAG_Node = it->second.CAG_Node_id;
        CagFind(CAG_Node, cag);
      }
      for (auto it = cag.cbegin(); it != cag.cend();) {
        CAG_Node CAG_Node_id = it->second.CAG_Node_id;
        CAG_Node current_root = it->second.pointing_to;

        auto find_it = roots_outgoing_edges_left.find(current_root);
        if (find_it == roots_outgoing_edges_left.end()) {
          if (DEBUG_CONDITION) {
            std::cout << "└─removing CAG_Node_id " << CAG_Node_id
                      << " because component has no outgoing edges. (root: "
                      << current_root << " )" << std::endl;
          }
          if (local_cc_size && CAG_Node_id >= lowest_local_CAG_Node_id &&
              CAG_Node_id <= highest_local_CAG_Node_id) {
            local_CAG_Nodes_pointer_table[CAG_Node_id].pointing_to =
                current_root;
          }
          it = cag.erase(it);
        } else {
          ++it;
        }
      }

      if (DEBUG_CONDITION) {
        std::cout << std::endl << "current cag:" << std::endl;
        for (auto it = cag.begin(); it != cag.end(); ++it) {
          std::cout << "└─" << it->second.ToString() << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
  if (sem("domain_overwrite")) {
    for (auto l : layers) {
      for (auto c : m.Cells()) {
        if (fccl_new[l][c] != kClNone) {
          (*fccl[l])[c] =
              1 +
              sin(local_CAG_Nodes_pointer_table[fccl_new[l][c]].pointing_to);
          //(*fccl[l])[c] =
          // local_CAG_Nodes_pointer_table[fccl_new[l][c]].pointing_to;
          // if (DEBUG_CONDITION2) std::cout << "no kCLNone: " << fccl_new[l][c]
          // << " " << local_CAG_Nodes_pointer_table[fccl_new[l][c]].pointing_to
          // << std::endl;
        } else {
          (*fccl[l])[c] = kClNone;
          // if (DEBUG_CONDITION2) std::cout << "nothin" << std::endl;
        }
      }
    }
  }
}