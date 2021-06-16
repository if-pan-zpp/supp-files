#pragma once
#include <string>
#include <vector>

namespace mdk::pdb {
    /**
     * An object representing a single field of a PDB record, along with
     * the means of parsing and printing it. The field objects signal parsing
     * failure by throwing.
     */
    class Field {
    public:
        /**
         * Read the value from a string (usually a line of a PDB file).
         * @param s String to read the value from.
         */
        virtual void read(const std::string &s) = 0;

        /**
         * Write the value to a string (usually a line of a PDB file).
         * Overflows are handled depending on the type of field.
         * @param s String to write the value to.
         */
        virtual void write(std::string &s) const = 0;
    };

    /**
     * Integer field.
     */
    class Integer: public Field {
    private:
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * End of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int j;

        /**
         * A value that is added to the read value and subtracted when a value
         * is written. It's chiefly used to convert 0-indexed ordinals to
         * 1-indexed and vice versa.
         */
        int offset;

        /**
         * A pointer to the integer to save the value to and load the value
         * from. It allows us in particular to have the same field class for
         * many different PDB records, only having to provide different pointer
         * to the underlying value.
         */
        int *v;

    public:
        Integer(int i, int j, int& v, int offset = 0):
            i{i}, j{j}, offset{offset}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    /**
     * Real (floating number) field.
     */
    class Real: public Field {
    private:
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * End of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int j;

        /**
         * Width of the read real number. Follows the conventions from
         * the PDB format docs and Fortran specifiers (it's n in f{n}.{m})
         */
        int n;

        /**
         * Precision of the read real number. Follows the conventions from
         * the PDB format docs and Fortran specifiers (it's m in f{n}.{m})
         */
        int m;

        /**
         * A pointer to the double to save the value to and load the value
         * from. It allows us in particular to have the same field class for
         * many different PDB records, only having to provide different pointer
         * to the underlying value.
         */
        double *v;

        /**
         * Scalar by which to multiply a read value and divide a written value.
         * In particular, it allows us to convert PDB file formats to internal
         * units and back.
         */
        double scalar;

    public:
        Real(int i, int j, int n, int m, double& v, double scalar = 1.0):
            i{i}, j{j}, n{n}, m{m}, v{&v}, scalar{scalar} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    /**
     * A mode enum, which modifies how to position and parse strings.
     * Note: most enums in the program are enum classes; here we need a standard
     * enum in order for combining via \p | to work
     */
    enum Mode {
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        Exact = 0x1,
        /**
         * When parsing, the final value is a trimmed string read from the
         * line.
         */
        Trim = 0x2,
        /**
         * When writing, the string is written left-aligned in the line.
         */
        Left = 0x4,
        /**
         * When writing, the string is written right-aligned in the line.
         */
        Right = 0x8
    };

    /**
     * String field. By default read values are trimmed, and written values
     * are right-aligned.
     */
    class String: public Field {
    private:
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * End of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int j;

        /**
         * A pointer to the string to save the value to and load the value
         * from. It allows us in particular to have the same field class for
         * many different PDB records, only having to provide different pointer
         * to the underlying value.
         */
        std::string *v;

        /**
         * Mode of string parsing and writing. See \p Mode for details.
         */
        int8_t mode;

    public:
        String(int i, int j, std::string& v, int8_t mode = Trim|Right):
            i{i}, j{j}, v{&v}, mode{mode} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    /**
     * (Single) character field.
     */
    class Char: public Field {
    private:
        /**
         * The location of the character in the line. Inclusive and 1-indexed,
         * as in the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * A pointer to the character to save the value to and load the value
         * from. It allows us in particular to have the same field class for
         * many different PDB records, only having to provide different pointer
         * to the underlying value.
         */
        char *v;

    public:
        Char(int i, char& v):
            i{i}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    /**
     * "Symmetry operator" field. At the moment it's parsed just like a normal
     * string.
     */
    class SymOp: public Field {
    private:
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * End of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int j;

        /**
         * A pointer to the string to save the value to and load the value
         * from. It allows us in particular to have the same field class for
         * many different PDB records, only having to provide different pointer
         * to the underlying value.
         */
        std::string *v;

    public:
        SymOp(int i, int j, std::string& v):
            i{i}, j{j}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    /**
     * A literal. As opposed to a normal field, here there's no reading or
     * writing values by pointers, rather comparing the read value with
     * the literal and printing out the literal. It's chiefly used for parsing
     * the record headers, like "ATOM  " or "MODEL ".
     */
    class Literal: public Field {
    private:
        /**
         * Beginning of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int i;

        /**
         * End of the field in the line. Inclusive and 1-indexed, as in
         * the PDB format docs, in order to make writing new fields easier.
         */
        int j;

        /**
         * Literal to compare the read value with and print out to the line
         * when writing.
         */
        std::string lit;

    public:
        Literal(int i, int j, std::string lit):
            i{i}, j{j}, lit{move(lit)} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };
}