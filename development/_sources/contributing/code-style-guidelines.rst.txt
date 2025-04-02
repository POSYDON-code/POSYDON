.. _code-style:

Code Style Guidelines
---------------------

Adhering to a consistent code style helps maintain readability and eases collaboration. Here are the guidelines you should follow when contributing to POSYDON.

General Principles
~~~~~~~~~~~~~~~~~~
1. **Clarity and Readability**: Code should be easy to understand. Opt for clarity over cleverness.
2. **Consistency**: Follow the style of the existing codebase. When in doubt, consistency with the existing code is more important than adhering to external style guides.
3. **Comments**: Write meaningful comments, especially for complex logic. Avoid obvious comments.

Python Style Guide
~~~~~~~~~~~~~~~~~~
We follow the `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide for Python code with some exceptions and additions as outlined below.

1. **Indentation**: Use 4 spaces per indentation level.
2. **Line Length**: Limit lines to 79 characters.
3. **Imports**: Group imports in the following order: standard libraries, third-party libraries, and local application imports. Each group should be separated by a blank line.
4. **Whitespace**: Avoid extraneous whitespace in the following situations:
    - Immediately inside parentheses, brackets, or braces.
    - Immediately before a comma, colon, or semicolon.
5. **Comments**: Comments should be complete sentences and should be used sparingly, i.e., only when necessary to explain complex pieces of logic or decisions that could seem non-obvious to other developers.

Naming Conventions
~~~~~~~~~~~~~~~~~~
1. **Functions and Variables**: Use `snake_case`.
2. **Classes and Exceptions**: Use `CamelCase`.
3. **Constants**: Use `UPPER_SNAKE_CASE`.
4. **Module-level globals** (like loggers or constants): Prefix with an underscore `_`.

Docstrings
~~~~~~~~~~
We adhere to the `NumPy docstring style guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_. Every function, class, and module should have a docstring that adheres to this style. This helps in generating consistent and readable documentation. Ensure:

1. Every module, class, and function/method has an accompanying docstring.
2. Always describe parameters, return types, and raised exceptions.
3. Include examples in function/method docstrings when possible.

Testing
~~~~~~~
Ensure that your code does not introduce new bugs:

1. Write unit tests for your new features.
2. Run all the tests to make sure they pass.
3. Follow the existing testing conventions in the codebase.

Commit Guidelines
~~~~~~~~~~~~~~~~~
1. Commit messages should be concise and descriptive.
2. Use the present tense ("Add feature" not "Added feature").
3. Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
4. Limit the first line to 72 characters or less.
5. Reference issues and pull requests liberally after the first line.

Final Thoughts
~~~~~~~~~~~~~~
Before submitting your code, review these guidelines and check your code against them. Taking the time to ensure your code adheres to these standards will make the review process smoother and faster.

Thank you for your contribution and for making POSYDON's codebase clean and consistent!
