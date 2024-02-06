# Pelops architecture

Documentation of the high level architecture design:
- module dependencies and relations
- class diagrams
- object diagrams
- factories

Files are created using [Umlet](https://www.umlet.com/), a free, open-source UML tool for which
there is also a Visual Studio Code extension and a web app called
[UMLetino](https://www.umletino.com/umletino.html) which can be used to load
and edit diagrams (`*.uxf` files).

The text version stores the latest correct design, but you can also see a png
version that can be visualised directly.

## Module dependencies and relations
This diagram provides the highest level view on all of the codebase. Shows how
different modules are related to each other, in particular how they depend on each other.

## Class diagrams
This diagrams shows the relationships of the most important classes in the
business logic: the callers and the components used by them. This helps to
provide a detailed understanding of the code by only understanding class
interfaces.

## Factories
The factories take care of building the relevant instance of objects depending
on the use case. They are not trivial as many classes have polymorphic
behaviours. This diagram will help to understand how objects are built.

## Object diagrams
This diagram shows how a few concrete objects are after construction.
The idea is not to show all possible different objects, but provide an
example for the default and probably most complex use case.
