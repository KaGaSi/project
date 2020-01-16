/**
 * \file
 * \brief Error prints
 */

#ifndef _ERRORS_H_
#define _ERRORS_H_

// ErrorCoorRead() //{{{
/** \brief Incorrect reading of vcf file
 *
 * \param [in] input_vcf  .vcf coordinate file
 * \param [in] bead       bead's line in its timestep in .vcf file where error occurred
 * \param [in] step       timestep when error occurred
 * \param [in] stuff      comment line of a timestep when error occurred
 * \param [in] input_vsf  .vsf structure file
 */
void ErrorCoorRead(char *input_vcf, int bead, int step, char *stuff, char *input_vsf); //}}}

// ErrorArgNumber() //{{{
/** \brief Insufficient number of arguments
 *
 * \param [in] count  number of supplied arguments
 * \param [in] need   minimum number of required arguments
 */
void ErrorArgNumber(int count, int need); //}}}

// ErrorOption() //{{{
/** \brief Unknown option
 *
 * \param [in] option   the unknown option
 */
void ErrorOption(char *option); //}}}

// ErrorNaN() //{{{
/** \brief Non-numeric argument
 *
 * \param [in] option   the option with wrong argument
 */
void ErrorNaN(char *option); //}}}

// ErrorExtension() //{{{
/** \brief Wrong file extension
 *
 * \param [in] file       filename
 * \param [in] number     number of correct extension(s)
 * \param [in] extension  correct extension(s)
 */
bool ErrorExtension(char *file, int number, char extension[][5]); //}}}

// ErrorFileOpen() //{{{
/** \brief Cannot open file
 *
 * \param [in] file  filename
 * \param [in] mode  open mode - r(ead), w(rite), a(ppend)
 */
void ErrorFileOpen(char *file, char mode); //}}}
#endif
