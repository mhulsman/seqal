#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "needleman_wunsch.h"


void nw_align(alignment_t * result, const char* seq_a, const char* seq_b, int match, int mismatch, int gap_open, int gap_extend, 
                            int no_start_gap_penalty, int no_end_gap_penalty, int no_gaps_in_a, int no_gaps_in_b, int no_mismatches, int case_sensitive)
{
  }


static PyObject * nw_align_wrapper(PyObject *self, PyObject *args, PyObject *kw)
{
    const char *seq1, *seq2;
    // Decide on scoring
    int match = 1;
    int mismatch = -2;
    int gap_open = -4;
    int gap_extend = -1;
    
    // Don't penalise gaps at the start
    // ACGATTT
    // ----TTT would score +3 (when match=+1)
    int no_start_gap_penalty = 0;
    
    // ..or gaps at the end e.g.
    // ACGATTT
    // ACGA--- would score +4 (when match=+1)
    int no_end_gap_penalty = 0;

    int no_gaps_in_a = 0, no_gaps_in_b = 0;
    int no_mismatches = 0;

    // Compare character case-sensitively (usually set to 0 for DNA etc)
    int case_sensitive = 0;

    PyObject * matrix = NULL;

    static char *kwlist[] = {"seq1","seq2", "matrix", "match", "mismatch", "gap_open","gap_extend", "no_start_gap_penalty", "no_end_gap_penalty", "no_gaps_in_a", "no_gaps_in_b", "no_mismatches", "case_sensitive", NULL};
    PyObject *res = NULL;

    if(!PyArg_ParseTupleAndKeywords(args, kw, "ss|Oiiiiiiiiii", kwlist, &seq1, &seq2, &matrix, &match, &mismatch, &gap_open, &gap_extend,
                                                                 &no_start_gap_penalty, &no_end_gap_penalty, &no_gaps_in_a, &no_gaps_in_b, &no_mismatches, &case_sensitive))
        return NULL;
    alignment_t *result = alignment_create(256);
    
    // Variables to store alignment result
    nw_aligner_t *nw = needleman_wunsch_new();

    scoring_t scoring;
    scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
                no_start_gap_penalty, no_end_gap_penalty,
                no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);

    // Add some special cases
    // x -> y means x in seq1 changing to y in seq2
    if(matrix != NULL)
    {
        PyObject * mapping = PyMapping_Items(matrix);
        if(mapping == NULL)
            goto error;
        int n = PySequence_Size(mapping);
        PyObject *item;
        int value;
        PyObject *key;
        char * char_a;
        char * char_b;
        int i;
        for(i = 0; i < n; i++)
        {
            item = PySequence_GetItem(mapping, i);
            if(item == NULL || !PyTuple_Check(item))
            {
                Py_XDECREF(item);
                Py_DECREF(mapping);
                goto error; 
            }
            
            if(!PyArg_ParseTuple(item, "Oi", &key, &value))
            {
                PyErr_SetString(PyExc_RuntimeError, "Values of matrix dict should be integers");
                Py_XDECREF(item);
                Py_DECREF(mapping);
                goto error;
            }
            if(!PyTuple_Check(key))
            {
                PyErr_SetString(PyExc_RuntimeError, "Keys of matrix dict should be tuples");
                Py_XDECREF(item);
                Py_DECREF(mapping);
                goto error;
            }
            if(!PyArg_ParseTuple(key, "ss", &char_a, &char_b))
            {
                PyErr_SetString(PyExc_RuntimeError, "Keys of matrix dict should be tuples with 2 characters as elements.");
                Py_XDECREF(item);
                Py_DECREF(mapping);
                goto error;
            }
            if(strlen(char_a) != 1 || strlen(char_b) != 1)
            {
                PyErr_SetString(PyExc_RuntimeError, "Character length should be 1");
                Py_XDECREF(item);
                Py_DECREF(mapping);
                goto error;
            }
            scoring_add_mutation(&scoring, case_sensitive ? *char_a : tolower(*char_a), case_sensitive ? *char_a : tolower(*char_b), value); // a -> c give substitution score -2
            Py_DECREF(item);
        }
    }

    // We could also prohibit the aligning of characters not given as special cases
    // scoring.use_match_mismatch = 0;

    needleman_wunsch_align(seq1, seq2, &scoring, nw, result);

    res = Py_BuildValue("ssi", result->result_a, result->result_b, result->score);

error:
    // Free memory for storing alignment results
    needleman_wunsch_free(nw);

    alignment_free(result);
    return res;
}

static PyMethodDef module_methods[] = {
    {"nw_align", (PyCFunction) nw_align_wrapper, METH_VARARGS | METH_KEYWORDS, "Needleman-Wunsch"},
    {NULL}  /* Sentinel */
};


#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initseqalign(void) 
{
    PyObject* m;
    m = Py_InitModule3("seqalign", module_methods,
                       "Sequence alignment utilities");
    return;

}


