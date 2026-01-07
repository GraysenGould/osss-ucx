/* For license: see LICENSE file at top-level */

/**
 * @file collectives.h
 * @brief Header file for OpenSHMEM collective operations
 *
 * This file contains declarations for initializing and finalizing
 * the collective operations subsystem.
*/

#ifndef _REDUCTIONS_H
#define _REDUCTIONS_H 1

/**
 * @brief Initialize the collective operations subsystem
 */
extern void collectives_init(void);

/**
 * @brief Finalize and cleanup the collective operations subsystem
 */
extern void collectives_finalize(void);

#endif /* ! _REDUCTIONS_H */
