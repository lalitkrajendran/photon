import numpy as np


def generate_progress_string(progress_decimal,progress_bar_character_number):
# % This function generates the progress bar text that will be displayed for
# % the current percentage of the calculation that is completed.

    # % This is a string giving a numerical value of the percentage of the
    # % calculation completed
    numerical_progress_string = '% 4.0f' % np.round(100*progress_decimal)

    # % This is the prefix to the progress bar showing the numerical value of the
    # % percent complete
    string_prefix = 'Calculation' + numerical_progress_string + '%% Complete:     0%% ['

    # % This is the suffix to the progress bar
    string_suffix = '] 100%%'

    # % This is the number of '#' signs to display corresponding to the graphical
    # % representation of the percent of the calculation that is complete
    completed_character_number = np.round(progress_decimal*progress_bar_character_number)

    # % This is the number of '-' signs to display corresponding to the graphical
    # % representation of the percent of the calculation that remains
    remaining_character_number = progress_bar_character_number - completed_character_number

    # % This creates the string of characters representing the graphical
    # % percentage of the calculation that is complete
    progress_bar_string = '#'*completed_character_number + '-'*remaining_character_number

    # % This is the total text string to be displayed
    text_string = string_prefix + progress_bar_string + string_suffix

    
    return text_string


def display_calculation_progress(current_value,value_vector):
# % This function displays a progress bar showing the percent complete of the
# % currently running calculation where the calculation is being iterated
# % over the vector 'value_vector' and the current iteration's value is equal
# % to 'current_value'.  For example in the for loop
# %
# %  for ii=1:10;
# %       commands . . .
# %  end;
# %
# % the arguments would be equal to
# %
# %  current_value=ii;
# %  value_vector=1:10;
# %
# % This convention holds for non-integer or non-monotonic vectors, although
# % to work correctly all values of value_vector must be unique.

    # % This is the number of characters to display in the pogress bar (denoted
    # % by either '#' or '-' characters)
    progress_bar_character_number = 50

    # % This is the index of the value_vector that corresponds to the
    # % current_value
    value_index = np.min(np.abs([current_value - x for x in value_vector]))

    # % This is the percentage of the calculation that is completed
    current_progress_decimal = value_index*1.0/len(value_vector)

    # % This creates the character string showing the numerical and graphical
    # % progress for the current iteration
    current_text_string = generate_progress_string(current_progress_decimal,progress_bar_character_number)

    # % If this is the first iteration, then a new line is added, the progress
    # % bar is displayed, and the function is exited
    if value_index==0:
        # % This displays the portion of the string showing the percentage of
        # % the calculation that is complete that is new to this iteration
        print current_text_string
        # % This ends the function
        return

    # % This is the percentage of the calculation that was completed during the
    # % last iteration
    previous_progress_decimal = (value_index - 1.0)/len(value_vector)

    # % This creates the character string showing the numerical and graphical
    # % progress for the previous iteration
    previous_text_string = generate_progress_string(previous_progress_decimal,progress_bar_character_number)

    # % This compares the current progress string with the previous string and if
    # % they are the same, then the function exits.  If they are different, then
    # % only the text after the difference is displayed.
    if current_text_string == previous_text_string:

        # % If this is the last time that the progress bar will be displayed, this
        # % prints a new line
        if value_index == length(value_vector)-1:
            # % This prints a new line after the progress bar
            print '\n'

        # % This exits the function without changing the progress bar
        return

    else:
        # % This is the total number of charcters to be displayed
        string_character_number = len(current_text_string)

        flag = 0
        # % This is the index into the string where the strings first differ
        for i in range(0,len(current_text_string)):
            if(current_text_string[:i+1] != previous_text_string[:i+1]):
                first_difference_index = i
                flag = 1

        if flag==0:
            first_difference_index = len(current_text_string)

        # % These are the locations of the double percent signs '%%'
        double_percent_indices = current_text_string.find('%%')

        # % This removes the double percent indices that are before the first
        # % difference index
        double_percent_indices = double_percent_indices[first_difference_index:]

        # % This is the number of characters of the previous line to delete
        delete_character_number = string_character_number - first_difference_index + 1 - len(double_percent_indices)

        # % If this is the first iteration, then no characters are deleted
        if value_index==1:
            # % This sets the number of characters to be deleted to zero
            delete_character_number = 0
            # % This sets the first difference character to one
            first_difference_index = 1

        # % This deletes the previously displayed characters back to the first
        # % differing character (by displaying the 'backspace' character)
        print '\b'*delete_character_number

        # % This displays the portion of the string showing the percentage of
        # % the calculation that is complete that is new to this iteration
        print current_text_string[first_difference_index:-1]

        # % If this is the last time that the progress bar will be displayed, this
        # % prints a new line
        if value_index==len(value_vector):
            # % This prints a new line after the progress bar
            print "\n"

